#!/usr/bin/env python3
"""
Stage 8 -- the systematic terms of every point, and save.

Reads 02_per_energy.csv (nominal rows) and, for the centroid pass, 07_uniformity.csv.
Two recipes, the two of the repository after merge #2; the selection decides which:

  codiceA (hodoscope, resolution_hodo.py)
      BES from colls_energies_summary_<R>ohm.csv: BES_cons is subtracted, BES_formula
      gives the "larger BES" variation
      sigma_corr      = sqrt(sigma^2 - BES_cons^2 - sync^2)
      bes_syst        = |sqrt(sigma^2 - BES_formula^2 - sync^2) - sigma_corr|
      syn_syst        = |sqrt(sigma^2 - BES_cons^2 - (1.3 sync)^2) - sigma_corr|
      err_total       = drift (+) stat (+) vtx_syst (+) bes_syst (+) syn_syst
  uniforme (centroid, resolution_final_uniforme.py --central raw --syst map)
      BES from rereco_<R>_withBES.csv (column bes); a point without BES is dropped
      sigma_corr      = sqrt(sigma^2 - BES^2 - sync^2)         (no position correction)
      unif_syst       = sigma * |s_energy - s_mean| / sigma_corr, the spread of the
                        corrected sigma/mu between the surface per energy and the surface
                        of the resistance, propagated through the subtraction
      err_total       = stat (+) drift (+) unif_syst
      pos_eff         = sqrt(sigma^2 - s_energy^2), written for information only

sync = 1.92e-7 * E_true^2.5 in percent. Every term is also written as a percentage of
the measured sigma/mu (the *_frac columns). Writes 08_systematics.csv.
"""

import argparse
import csv
import math
import os

import numpy as np

import common

COLUMNS = ("resistance", "energy", "energy_true", "selection", "recipe", "window", "n_events", "n_run",
           "pooled", "sigma_raw", "stat", "drift", "chi2_drift", "syst_tails", "vtx_syst", "bes", "bes_nom",
           "sync", "bes_syst", "syn_syst", "unif_syst", "pos_eff", "pos_naive", "s_run", "s_energy",
           "s_mean", "err_total", "sigma_corr",
           "stat_frac", "drift_frac", "tails_frac", "vtx_frac", "bes_frac", "sync_frac", "unif_frac")


def load_bes_codice_a(besdir, resistance):
    """{energy: (BES_cons, BES_formula)} from colls_energies_summary_<R>ohm.csv."""
    path = os.path.join(besdir, f"colls_energies_summary_{resistance}ohm.csv")
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found: BES = 0 for {resistance} ohm (resolution_hodo.py does the same)")
        return {}
    out = {}
    with open(path) as handle:
        for row in csv.DictReader(handle):
            out[int(float(row["Energy"]))] = (float(row["BES_cons"]), float(row["BES_formula"]))
    return out


def load_bes_uniforme(besdir, resistance):
    """{energy: bes} from rereco_<R>_withBES.csv, column 7 (uniformita_pos.load_bes)."""
    path = os.path.join(besdir, f"rereco_{resistance}_withBES.csv")
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found: every {resistance} ohm point is dropped (no BES)")
        return {}
    out = {}
    with open(path) as handle:
        for index, line in enumerate(handle):
            if index == 0:
                continue
            parts = line.strip().split(",")
            if len(parts) >= 7:
                out[int(float(parts[0]))] = float(parts[6])
    return out


def subtract(sigma, *terms):
    return math.sqrt(max(sigma * sigma - sum(term * term for term in terms), 0.))


def codice_a_row(point, bes_table):
    sigma = point["sigma_over_mu"]
    bes_cons, bes_formula = bes_table.get(point["energy"], (0., 0.))
    sync = common.synchrotron_pct(point["energy_true"])
    sigma_corr = subtract(sigma, bes_cons, sync)
    bes_syst = abs(subtract(sigma, bes_formula, sync) - sigma_corr)
    syn_syst = abs(subtract(sigma, bes_cons, 1.3 * sync) - sigma_corr)
    err_total = math.sqrt(point["drift"] ** 2 + point["stat"] ** 2 + point["vtx_syst"] ** 2
                          + bes_syst ** 2 + syn_syst ** 2)
    return dict(bes=bes_cons, bes_nom=bes_formula, sync=sync, bes_syst=bes_syst, syn_syst=syn_syst,
                unif_syst=np.nan, pos_eff=np.nan, pos_naive=np.nan, s_run=np.nan, s_energy=np.nan,
                s_mean=np.nan, err_total=err_total, sigma_corr=sigma_corr if sigma_corr > 0 else np.nan)


def uniforme_row(point, bes_table, uniformity):
    sigma = point["sigma_over_mu"]
    bes = bes_table.get(point["energy"], 0.)
    if not bes > 0:
        print(f"  {point['resistance']} ohm {point['energy']:>4} GeV: dropped, no BES in rereco_{point['resistance']}_withBES.csv")
        return None
    sync = common.synchrotron_pct(point["energy_true"])
    sigma_corr = subtract(sigma, bes, sync)
    map_syst = uniformity.get("syst_pct", 0.) if uniformity else 0.
    if not np.isfinite(map_syst):
        map_syst = 0.
    unif_syst = sigma * map_syst / sigma_corr if sigma_corr > 0 else 0.
    err_total = math.sqrt(point["stat"] ** 2 + point["drift"] ** 2 + unif_syst ** 2)
    s_energy = uniformity.get("s_energy", np.nan) if uniformity else np.nan
    return dict(bes=bes, bes_nom=bes, sync=sync, bes_syst=np.nan, syn_syst=np.nan, unif_syst=unif_syst,
                pos_eff=subtract(sigma, s_energy) if np.isfinite(s_energy) else np.nan,
                pos_naive=uniformity.get("pos_term", np.nan) if uniformity else np.nan,
                s_run=uniformity.get("s_run", np.nan) if uniformity else np.nan, s_energy=s_energy,
                s_mean=uniformity.get("s_mean", np.nan) if uniformity else np.nan,
                err_total=err_total, sigma_corr=sigma_corr if sigma_corr > 0 else np.nan)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--besdir", required=True, help="directory with the BES tables")
    args = parser.parse_args()

    points = [row for row in common.read_csv(common.require(os.path.join(args.workdir, "02_per_energy.csv")))
              if row["variation"] == "nominal"]
    if not points:
        raise SystemExit("no nominal point in 02_per_energy.csv")
    recipe = points[0]["recipe"]
    uniformity = {}
    if recipe == "uniforme":
        uniformity = {(row["resistance"], row["energy"]): row
                      for row in common.read_csv(common.require(os.path.join(args.workdir, "07_uniformity.csv")))}

    rows = []
    bes_tables = {}
    for point in sorted(points, key=lambda row: (row["resistance"], row["energy"])):
        resistance = point["resistance"]
        if resistance not in bes_tables:
            bes_tables[resistance] = (load_bes_codice_a if recipe == "codiceA" else load_bes_uniforme)(args.besdir, resistance)
        if recipe == "codiceA":
            terms = codice_a_row(point, bes_tables[resistance])
        else:
            terms = uniforme_row(point, bes_tables[resistance], uniformity.get((resistance, point["energy"])))
        if terms is None:
            continue
        row = dict(resistance=resistance, energy=point["energy"], energy_true=point["energy_true"],
                   selection=point["selection"], recipe=recipe, window=point["window"],
                   n_events=point["n_events"], n_run=point["n_run"], pooled=point["pooled"],
                   sigma_raw=point["sigma_over_mu"], stat=point["stat"], drift=point["drift"],
                   chi2_drift=point["chi2_ndf_drift"], syst_tails=point["syst_tails"], vtx_syst=point["vtx_syst"])
        row.update(terms)
        sigma = row["sigma_raw"]
        for name, source in (("stat_frac", "stat"), ("drift_frac", "drift"), ("tails_frac", "syst_tails"),
                             ("vtx_frac", "vtx_syst"), ("bes_frac", "bes"), ("sync_frac", "sync"),
                             ("unif_frac", "unif_syst")):
            row[name] = 100. * row[source] / sigma if np.isfinite(row[source]) and sigma > 0 else np.nan
        rows.append(row)
        print(f"  {resistance} ohm {point['energy']:>4} GeV: sigma/mu {sigma:.4f}  -> sigma/E {row['sigma_corr']:.4f} "
              f"+- {row['err_total']:.4f} %   (stat {row['stat']:.4f}, drift {row['drift']:.4f}, BES {row['bes']:.3f}, "
              f"sync {row['sync']:.3f}"
              + (f", vtx {row['vtx_syst']:.4f}, BES syst {row['bes_syst']:.4f}, sync syst {row['syn_syst']:.4f})"
                 if recipe == "codiceA" else f", map syst {row['unif_syst']:.4f})"))
    common.write_csv(os.path.join(args.workdir, "08_systematics.csv"), rows, COLUMNS)


if __name__ == "__main__":
    main()
