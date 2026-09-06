#!/usr/bin/env python3
"""
Stage 6 -- correct the amplitude event by event for the response surface, and save.

Run once per --mode, which chooses the surface f(u, v) used for each run
(uniformita_maps.py, stage apply):

  run      the surface of the run itself (04_surface_per_run.csv); a run without one
           (fewer than 200 events in its window) falls back on the surface of its energy
  energy   the surface of the energy, all its runs together (05_surface_mean.csv)
  mean     the single surface of the resistance, all energies together

The correction is A -> A * <f>_window / f(u_i, v_i), with <f>_window the mean of f over
the events of the run inside its fit window, so the peak of the run does not move. The
position term of uniformita_pos.py, 100 * std(f) / <f> over the same events, is saved
next to it (pos_pct).

Writes, in --workdir:
  06_corrected_<mode>.root   TTree "corrected" with resistance, energy, run, amplitude and
                             corrected amplitude per event (stage 7 fits it with the same
                             window and binning recipe as stage 1), plus one TH1D per run
                             and per energy of the corrected amplitude, 1 ADC bins
  06_correction_<mode>.csv   per run: which surface was used, <f>, pos_pct
"""

import argparse
import os

import numpy as np
import ROOT

import common
from common import runsets
from s03_maps_and_profiles import design

COLUMNS = ("resistance", "energy", "run", "mode", "surface", "fallback", "n_selected", "n_window",
           "f_mean", "pos_pct")


def coefficients_of(row):
    return np.array([row[f"a{index}"] for index in range(6)], float)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--mode", choices=("run", "energy", "mean"), required=True)
    parser.add_argument("--resistances", nargs="+", type=int, default=[340, 400, 500])
    parser.add_argument("--exclude-runs", nargs="*", type=int, default=[])
    runsets.add_argument(parser)
    args = parser.parse_args()

    fit_rows = common.read_csv(common.require(os.path.join(args.workdir, "01_dcb_per_run.csv")))
    windows = {(row["resistance"], row["energy"]): row
               for row in common.read_csv(common.require(os.path.join(args.workdir, "01_windows.csv")))}
    moment_runs = {(row["resistance"], row["energy"], row["run"]): row
                   for row in common.read_csv(common.require(os.path.join(args.workdir, "03_moments.csv")))}
    run_surfaces = {(row["resistance"], row["energy"], row["run"]): row
                    for row in common.read_csv(common.require(os.path.join(args.workdir, "04_surface_per_run.csv")))}
    mean_surfaces = {}
    for row in common.read_csv(common.require(os.path.join(args.workdir, "05_surface_mean.csv"))):
        mean_surfaces[(row["resistance"], row["energy"], row["scope"])] = row
    amplitude_name = fit_rows[0]["amplitude"] if fit_rows else "a3x3"
    dropped, kept_only = runsets.resolve(args.runset, args.exclude_runs)
    excluded = {key for key in windows if windows[key]["skipped"]}

    output = ROOT.TFile(os.path.join(args.workdir, f"06_corrected_{args.mode}.root"), "RECREATE")
    columns = {"resistance": [], "energy": [], "run": [], "amplitude": [], "corrected": []}
    rows = []
    for resistance, energy, path in common.resistance_energy_pairs(args.base, args.resistances, excluded):
        if (resistance, energy) not in windows:
            continue
        fits = {row["run"]: row for row in fit_rows
                if row["resistance"] == resistance and row["energy"] == energy
                and row["variation"] == "nominal" and row["tails_mode"] == "free" and row["run"] > 0
                and (resistance, energy, row["run"]) in moment_runs}
        energy_surface = mean_surfaces.get((resistance, energy, "energy"))
        resistance_surface = mean_surfaces.get((resistance, 0, "mean"))
        if not fits or energy_surface is None or resistance_surface is None:
            print(f"[{resistance} ohm {energy:>4} GeV] no surface, skipped")
            continue
        print(f"[{resistance} ohm {energy:>4} GeV] {len(fits)} runs, mode {args.mode}", flush=True)
        events = common.read_events(path, amplitude_name)
        half = windows[(resistance, energy)]["half"]
        selected = ((events["A_tot"] > common.A_TOT_MIN) & common.runset_mask(events["run"], dropped, kept_only)
                    & (np.abs(events["u"]) <= half) & (np.abs(events["v"]) <= half))
        energy_histogram = common.fill_histogram(f"h_{resistance}_{energy}_all", [], common.HISTOGRAM_NBINS,
                                                 common.HISTOGRAM_LO, common.HISTOGRAM_HI)
        for this_run, fit in sorted(fits.items()):
            in_run = selected & (events["run"] == this_run)
            if not in_run.any():
                continue
            surface, fallback = args.mode, 0
            if args.mode == "run":
                own = run_surfaces.get((resistance, energy, this_run))
                if own is not None and own["ok"]:
                    coefficients = coefficients_of(own)
                else:
                    coefficients, surface, fallback = coefficients_of(energy_surface), "energy", 1
            elif args.mode == "energy":
                coefficients = coefficients_of(energy_surface)
            else:
                coefficients = coefficients_of(resistance_surface)
            amplitude = events["amplitude"][in_run]
            response = design(events["u"][in_run], events["v"][in_run]) @ coefficients
            in_window = (amplitude >= fit["window_lo"]) & (amplitude <= fit["window_hi"])
            response_mean = float(response[in_window].mean()) if in_window.any() else np.nan
            if not (response_mean > 0):
                print(f"    run {this_run}: <f> not positive, skipped")
                continue
            corrected = amplitude * response_mean / np.maximum(response, 1e-9)
            rows.append(dict(resistance=resistance, energy=energy, run=this_run, mode=args.mode,
                             surface=surface, fallback=fallback, n_selected=int(in_run.sum()),
                             n_window=int(in_window.sum()), f_mean=response_mean,
                             pos_pct=100 * float(response[in_window].std()) / response_mean))
            columns["resistance"].append(np.full(len(amplitude), resistance, np.int32))
            columns["energy"].append(np.full(len(amplitude), energy, np.int32))
            columns["run"].append(np.full(len(amplitude), this_run, np.int32))
            columns["amplitude"].append(amplitude)
            columns["corrected"].append(corrected)
            histogram = common.fill_histogram(f"h_{resistance}_{energy}_{this_run}", corrected,
                                              common.HISTOGRAM_NBINS, common.HISTOGRAM_LO, common.HISTOGRAM_HI)
            histogram.SetTitle(f"{resistance} ohm {energy} GeV run {this_run}, corrected ({surface} surface)")
            output.cd()
            histogram.Write()
            energy_histogram.Add(histogram)
        energy_histogram.SetTitle(f"{resistance} ohm {energy} GeV all runs, corrected")
        output.cd()
        energy_histogram.Write()
    output.Close()

    if columns["run"]:
        arrays = {name: np.concatenate(values) for name, values in columns.items()}
        frame = ROOT.RDF.FromNumpy(arrays)
        options = ROOT.RDF.RSnapshotOptions()
        options.fMode = "UPDATE"
        frame.Snapshot("corrected", os.path.join(args.workdir, f"06_corrected_{args.mode}.root"),
                       list(arrays), options)
        print("->", os.path.join(args.workdir, f"06_corrected_{args.mode}.root"),
              f"({len(arrays['run'])} events)")
    common.write_csv(os.path.join(args.workdir, f"06_correction_{args.mode}.csv"), rows, COLUMNS)


if __name__ == "__main__":
    main()
