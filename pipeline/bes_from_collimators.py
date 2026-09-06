#!/usr/bin/env python3
"""
The beam energy spread table of resolution_hodo.py, colls_energies_summary_<R>ohm.csv,
rebuilt from the collimator log.

BES (RMS, in percent) from the half-widths of the two momentum-defining collimators
of H4, C3 and C8 (CERN-SL-Note-97-81, eq. 2, and the slides of 10 September 2026):

    BES_formula = sqrt(C3^2 + C8^2) / (27 sqrt(3))       the nominal formula
    BES_cons    = C3 / (27 sqrt(3))                       the lower bound, "conservative"

The formula holds for C8 << C3; when C8 ~ C3 the simulation deviates from it, so the
lower bound is what gets subtracted and the difference between the two is carried as a
systematic (see stage 8).

Inputs:
  --collimators   the xlsx export of the collimator jaws against UTC time
                  (CMS_ECAL_Collimators_June26.xlsx). C3 is XCHV.022.131, C8 is
                  XCSV.022.386, the collimators ordered by position along the line.
  --timestamps    timestamps_runs.txt, lines "Jun 12 00:37 20769" in CERN local time
                  (CEST = UTC + 2), the time stamp of the run file
  --runs-csv      any 01_dcb_per_run.csv: which runs belong to which (resistance, energy)

For every run the jaws are read at the run time stamp and --lookback minutes earlier; a
run that sits across a collimator change is flagged. The table per (resistance, energy)
takes the mean over its runs and flags the energies whose runs had different settings.

Writes, in --outdir: colls_energies_summary_<R>ohm.csv (Energy, C3, C8, BES_formula,
BES_cons, n_runs, consistent) and colls_per_run.csv.
"""

import argparse
import datetime
import math
import os
import re
import xml.etree.ElementTree as ElementTree
import zipfile
from collections import defaultdict

import numpy as np

import common

XML_NAMESPACE = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
C3_COLUMN_NAME, C8_COLUMN_NAME = "XCHV.022.131", "XCSV.022.386"
LOCAL_TIME_OFFSET = datetime.timedelta(hours=2)          # CEST
BES_DENOMINATOR = 27. * math.sqrt(3.)
YEAR = 2026
SUMMARY_COLUMNS = ("Energy", "C3", "C8", "BES_formula", "BES_cons", "n_runs", "consistent")
RUN_COLUMNS = ("resistance", "energy", "run", "time_utc", "C3", "C8", "C3_before", "C8_before",
               "stable", "BES_formula", "BES_cons")


def read_collimator_log(path):
    """[(utc datetime, C3 half-width, C8 half-width)] sorted in time."""
    archive = zipfile.ZipFile(path)
    sheet = ElementTree.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    shared = []
    if "xl/sharedStrings.xml" in archive.namelist():
        strings = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
        for item in strings.findall("m:si", XML_NAMESPACE):
            shared.append("".join(text.text or "" for text in item.iter("{%s}t" % XML_NAMESPACE["m"])))

    def cell_values(row):
        values = []
        for cell in row:
            value = cell.find("m:v", XML_NAMESPACE)
            text = value.text if value is not None else ""
            if cell.get("t") == "s" and text:
                text = shared[int(text)]
            values.append(text.strip())
        return values

    rows = list(sheet.find("m:sheetData", XML_NAMESPACE))
    header = cell_values(rows[0])
    c3_jaw2 = next(index for index, name in enumerate(header) if C3_COLUMN_NAME in name and "JAW2" in name)
    c8_jaw2 = next(index for index, name in enumerate(header) if C8_COLUMN_NAME in name and "JAW2" in name)
    log = []
    for row in rows[1:]:
        values = cell_values(row)
        if not values or not values[0]:
            continue
        when = datetime.datetime(1899, 12, 30) + datetime.timedelta(days=float(values[0]))
        log.append((when, float(values[c3_jaw2]), float(values[c8_jaw2])))
    return sorted(log)


def read_run_times(path):
    """{run: utc datetime} from lines like 'Jun 12 00:37 20769'."""
    out = {}
    for line in open(path):
        parts = line.split()
        if len(parts) >= 4 and parts[-1].isdigit():
            local = datetime.datetime.strptime(f"{YEAR} {' '.join(parts[:3])}", "%Y %b %d %H:%M")
            out[int(parts[-1])] = local - LOCAL_TIME_OFFSET
    return out


def setting_at(log, when):
    """The last jaw setting logged before `when`."""
    current = None
    for entry in log:
        if entry[0] <= when:
            current = entry
        else:
            break
    return current


def bes_values(c3_half, c8_half):
    return (math.sqrt(c3_half ** 2 + c8_half ** 2) / BES_DENOMINATOR, c3_half / BES_DENOMINATOR)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--collimators", required=True, help="xlsx export of the collimator jaws")
    parser.add_argument("--timestamps", required=True, help="timestamps_runs.txt")
    parser.add_argument("--runs-csv", required=True, nargs="+",
                        help="one or more 01_dcb_per_run.csv, for the run -> (resistance, energy) map")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--lookback", type=float, default=20., help="minutes before the time stamp")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    log = read_collimator_log(args.collimators)
    run_times = read_run_times(args.timestamps)
    runs = {}
    for path in args.runs_csv:
        for row in common.read_csv(path):
            if row["run"] > 0:
                runs[row["run"]] = (row["resistance"], row["energy"])

    run_rows = []
    for run, (resistance, energy) in sorted(runs.items()):
        if run not in run_times:
            print(f"  run {run}: no time stamp, skipped")
            continue
        when = run_times[run]
        now, before = setting_at(log, when), setting_at(log, when - datetime.timedelta(minutes=args.lookback))
        if now is None:
            print(f"  run {run}: before the first collimator entry, skipped")
            continue
        formula, conservative = bes_values(now[1], now[2])
        stable = int(before is not None and before[1:] == now[1:])
        if not stable:
            print(f"  run {run} ({resistance} ohm {energy} GeV): jaws changed within {args.lookback:.0f} min "
                  f"before {when:%b %d %H:%M} UTC: C3/C8 {before[1:] if before else None} -> {now[1:]}")
        run_rows.append(dict(resistance=resistance, energy=energy, run=run, time_utc=f"{when:%Y-%m-%d %H:%M}",
                             C3=now[1], C8=now[2], C3_before=before[1] if before else np.nan,
                             C8_before=before[2] if before else np.nan, stable=stable,
                             BES_formula=formula, BES_cons=conservative))
    common.write_csv(os.path.join(args.outdir, "colls_per_run.csv"), run_rows, RUN_COLUMNS)

    grouped = defaultdict(list)
    for row in run_rows:
        grouped[(row["resistance"], row["energy"])].append(row)
    for resistance in sorted({key[0] for key in grouped}):
        summary = []
        for (this_resistance, energy), rows in sorted(grouped.items()):
            if this_resistance != resistance:
                continue
            c3_values, c8_values = [row["C3"] for row in rows], [row["C8"] for row in rows]
            consistent = int(len(set(c3_values)) == 1 and len(set(c8_values)) == 1)
            c3_mean, c8_mean = float(np.mean(c3_values)), float(np.mean(c8_values))
            formula, conservative = bes_values(c3_mean, c8_mean)
            summary.append(dict(Energy=energy, C3=c3_mean, C8=c8_mean, BES_formula=formula,
                                BES_cons=conservative, n_runs=len(rows), consistent=consistent))
            note = "" if consistent else f"   <-- runs differ: C3 {sorted(set(c3_values))}, C8 {sorted(set(c8_values))}"
            print(f"  {resistance} ohm {energy:>4} GeV: C3 {c3_mean:4.1f}  C8 {c8_mean:4.1f} mm  "
                  f"BES formula {formula:.4f}  cons {conservative:.4f} %  ({len(rows)} runs){note}")
        common.write_csv(os.path.join(args.outdir, f"colls_energies_summary_{resistance}ohm.csv"), summary,
                         SUMMARY_COLUMNS)


if __name__ == "__main__":
    main()
