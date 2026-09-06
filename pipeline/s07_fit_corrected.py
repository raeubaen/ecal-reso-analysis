#!/usr/bin/env python3
"""
Stage 7 -- double Crystal Ball fit of the corrected amplitudes, and save.

For every 06_corrected_<mode>.root found (modes run, energy, mean) the corrected
amplitude of each run is fitted with the same recipe as stage 1 (fit.sh window,
Freedman-Diaconis binning, free tails, MIGRAD + HESSE), and the runs of each energy are
combined with weights = number of selected events, as uniformita_maps.py does.

Writes:
  07_dcb_corrected.csv        one row per (run, mode)
  07_dcb_corrected_fits.root  histograms and TF1 of every fit
  07_uniformity.csv           one row per (resistance, energy), the columns of
                              uniformita_maps.csv: raw (from stage 2), s_run, s_energy,
                              s_mean with their errors, the scatter of the run variant,
                              syst_pct = |s_energy - s_mean| (the centroid systematic of
                              the post-merge code), the curvatures, and pos_term (the
                              position term of uniformita_pos.py, energy surface)
  uniformity_<R>ohm.png       raw and the three corrected sigma/mu against the energy
"""

import argparse
import glob
import os
import re
from collections import defaultdict

import numpy as np
import ROOT

import common

MODES = ("run", "energy", "mean")
MODE_COLOUR = {"run": ROOT.kAzure + 2, "energy": ROOT.kRed + 1, "mean": ROOT.kGreen + 2}
MODE_MARKER = {"run": 21, "energy": 22, "mean": 23}
FIT_COLUMNS = ("resistance", "energy", "run", "mode", "surface", "fallback", "n_selected", "n_events",
               "peak", "err_peak", "sigma", "err_sigma", "sigma_over_mu", "err_sigma_over_mu",
               "chi2", "ndf", "alpha_l", "alpha_h", "n_l", "n_h")
UNIFORMITY_COLUMNS = (["resistance", "energy", "energy_true", "nrun", "nev", "n_fallback", "raw", "err_raw",
                       "scat_run"] + [f"{prefix}_{mode}" for mode in MODES for prefix in ("s", "err")]
                      + ["syst_pct", "syst_rel_pct", "pos_term", "curv_eta_energy", "curv_phi_energy",
                         "curv_eta_mean", "curv_phi_mean"])


def fit_mode(path, mode, correction_rows, output):
    frame = ROOT.RDataFrame("corrected", path)
    arrays = frame.AsNumpy(["resistance", "energy", "run", "corrected"])
    resistance_of, energy_of, run_of = (np.asarray(arrays[name]) for name in ("resistance", "energy", "run"))
    corrected = np.asarray(arrays["corrected"], float)
    rows = []
    for correction in correction_rows:
        resistance, energy, this_run = correction["resistance"], correction["energy"], correction["run"]
        in_run = (resistance_of == resistance) & (energy_of == energy) & (run_of == this_run)
        fit = common.fit_dcb(corrected[in_run], energy, resistance)
        if fit is None:
            print(f"    {resistance} ohm {energy} GeV run {this_run} ({mode}): fit failed")
            continue
        value, error = common.relative_width(fit)
        row = dict(resistance=resistance, energy=energy, run=this_run, mode=mode, surface=correction["surface"],
                   fallback=correction["fallback"], n_selected=correction["n_selected"],
                   n_events=fit["n_events"], peak=fit["peak"], err_peak=fit["err_peak"], sigma=fit["sigma"],
                   err_sigma=fit["err_sigma"], sigma_over_mu=value, err_sigma_over_mu=error,
                   chi2=fit["chi2"], ndf=fit["ndf"])
        row.update(fit["tails"])
        rows.append(row)
        output.cd()
        fit["histogram"].SetName(f"h_{resistance}_{energy}_{this_run}_{mode}")
        fit["histogram"].Write()
        fit["function"].SetName(f"f_{resistance}_{energy}_{this_run}_{mode}")
        fit["function"].Write()
    return rows


def draw_uniformity(rows, resistance, outdir):
    rows = sorted([row for row in rows if row["resistance"] == resistance], key=lambda row: row["energy"])
    if len(rows) < 2:
        return
    energies = np.array([row["energy_true"] for row in rows], float)
    canvas = ROOT.TCanvas(f"uniformity_{resistance}", "", 1000, 1200)
    canvas.Divide(1, 3)
    pad = canvas.cd(1)
    pad.SetLogx()
    pad.SetGrid()
    raw = common.keep(common.make_graph(energies, [row["raw"] for row in rows]))
    raw.SetMarkerStyle(20)
    raw.SetMarkerColor(ROOT.kGray + 2)
    raw.SetLineColor(ROOT.kGray + 2)
    raw.SetTitle(f"{resistance} #Omega, A_{{tot}} > {common.A_TOT_MIN:.0f} ADC, |pos_eta - 18| #leq 0.2, "
                 f"|pos_phi - 6| #leq 0.2;E_{{true}} [GeV];#sigma/#mu  [%]")
    raw.Draw("APL")
    raw.GetXaxis().SetLimits(15, 320)
    all_values = [row["raw"] for row in rows] + [row[f"s_{mode}"] for mode in MODES for row in rows]
    raw.GetYaxis().SetRangeUser(0.8 * np.nanmin(all_values), 1.15 * np.nanmax(all_values))
    legend = common.keep(ROOT.TLegend(0.6, 0.65, 0.88, 0.88))
    legend.AddEntry(raw, "no correction", "pl")
    for mode in MODES:
        graph = common.keep(common.make_graph(energies, [row[f"s_{mode}"] for row in rows],
                                              [row[f"err_{mode}"] for row in rows]))
        graph.SetMarkerStyle(MODE_MARKER[mode])
        graph.SetMarkerColor(MODE_COLOUR[mode])
        graph.SetLineColor(MODE_COLOUR[mode])
        graph.Draw("PL")
        legend.AddEntry(graph, f"surface per {mode}", "pl")
    legend.Draw()

    pad = canvas.cd(2)
    pad.SetLogx()
    pad.SetGrid()
    reference = np.array([row["s_energy"] for row in rows])
    ratios = [100 * (np.array([row[f"s_{mode}"] for row in rows]) / reference - 1) for mode in MODES]
    span = 1.3 * max(np.nanmax(np.abs(ratio)) for ratio in ratios) + 0.1
    first = True
    for mode, ratio in zip(MODES, ratios):
        graph = common.keep(common.make_graph(energies, ratio))
        graph.SetMarkerStyle(MODE_MARKER[mode])
        graph.SetMarkerColor(MODE_COLOUR[mode])
        graph.SetLineColor(MODE_COLOUR[mode])
        graph.SetTitle(";E_{true} [GeV];relative to the surface per energy  [%]")
        graph.Draw("APL" if first else "PL")
        if first:
            graph.GetXaxis().SetLimits(15, 320)
            graph.GetYaxis().SetRangeUser(-span, span)
        first = False

    pad = canvas.cd(3)
    pad.SetLogx()
    pad.SetGrid()
    syst = common.keep(common.make_graph(energies, [row["syst_pct"] for row in rows]))
    syst.SetMarkerStyle(33)
    syst.SetMarkerColor(ROOT.kViolet)
    syst.SetLineColor(ROOT.kViolet)
    syst.SetTitle(";E_{true} [GeV];|per energy - mean|  [%]")
    syst.Draw("APL")
    syst.GetXaxis().SetLimits(15, 320)
    common.save_canvas(canvas, os.path.join(outdir, f"uniformity_{resistance}ohm.png"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()
    common.style()

    per_energy = {(row["resistance"], row["energy"]): row
                  for row in common.read_csv(common.require(os.path.join(args.workdir, "02_per_energy.csv")))
                  if row["variation"] == "nominal"}
    surfaces = {(row["resistance"], row["energy"], row["scope"]): row
                for row in common.read_csv(common.require(os.path.join(args.workdir, "05_surface_mean.csv")))}

    output = ROOT.TFile(os.path.join(args.workdir, "07_dcb_corrected_fits.root"), "RECREATE")
    fit_rows = []
    for path in sorted(glob.glob(os.path.join(args.workdir, "06_corrected_*.root"))):
        mode = re.search(r"06_corrected_(\w+)\.root", path).group(1)
        correction_rows = common.read_csv(os.path.join(args.workdir, f"06_correction_{mode}.csv"))
        print(f"[{mode}] {len(correction_rows)} runs", flush=True)
        fit_rows += fit_mode(path, mode, correction_rows, output)
    output.Close()
    common.write_csv(os.path.join(args.workdir, "07_dcb_corrected.csv"), fit_rows, FIT_COLUMNS)

    grouped = defaultdict(list)
    for row in fit_rows:
        grouped[(row["resistance"], row["energy"])].append(row)
    uniformity_rows = []
    for (resistance, energy), rows in sorted(grouped.items()):
        raw = per_energy.get((resistance, energy))
        if raw is None:
            continue
        row = dict(resistance=resistance, energy=energy, energy_true=common.true_energy(energy),
                   raw=raw["sigma_over_mu"], err_raw=raw["err_mean"], nrun=raw["n_run"], nev=raw["n_events"],
                   n_fallback=sum(item["fallback"] for item in rows if item["mode"] == "run"))
        for mode in MODES:
            mode_rows = [item for item in rows if item["mode"] == mode]
            values = [item["sigma_over_mu"] for item in mode_rows]
            errors = [item["err_sigma_over_mu"] for item in mode_rows]
            weights = [item["n_selected"] for item in mode_rows]
            row[f"s_{mode}"], row[f"err_{mode}"] = common.weighted_mean(values, errors, weights)
            if mode == "run":
                row["scat_run"] = common.weighted_scatter(values, weights)
        row["syst_pct"] = (abs(row["s_energy"] - row["s_mean"])
                           if np.isfinite(row["s_energy"]) and np.isfinite(row["s_mean"]) else np.nan)
        row["syst_rel_pct"] = 100 * row["syst_pct"] / row["raw"] if row["raw"] > 0 else np.nan
        energy_rows = [item for item in rows if item["mode"] == "energy"]
        position_terms = common.read_csv(os.path.join(args.workdir, "06_correction_energy.csv")) \
            if os.path.exists(os.path.join(args.workdir, "06_correction_energy.csv")) else []
        position_terms = [item for item in position_terms
                          if item["resistance"] == resistance and item["energy"] == energy]
        row["pos_term"] = (common.weighted_mean([item["pos_pct"] for item in position_terms],
                                                [0.] * len(position_terms),
                                                [item["n_selected"] for item in position_terms])[0]
                           if position_terms else np.nan)
        energy_surface = surfaces.get((resistance, energy, "energy"), {})
        mean_surface = surfaces.get((resistance, 0, "mean"), {})
        row.update(curv_eta_energy=energy_surface.get("curv_eta_pct", np.nan),
                   curv_phi_energy=energy_surface.get("curv_phi_pct", np.nan),
                   curv_eta_mean=mean_surface.get("curv_eta_pct", np.nan),
                   curv_phi_mean=mean_surface.get("curv_phi_pct", np.nan))
        uniformity_rows.append(row)
        print(f"  {resistance} ohm {energy:>4} GeV: raw {row['raw']:.4f}  run {row['s_run']:.4f}  "
              f"energy {row['s_energy']:.4f}  mean {row['s_mean']:.4f}  |energy - mean| {row['syst_pct']:.4f}"
              + (f"  ({row['n_fallback']} run without own surface)" if row["n_fallback"] else ""))
    common.write_csv(os.path.join(args.workdir, "07_uniformity.csv"), uniformity_rows, UNIFORMITY_COLUMNS)
    for resistance in sorted({row["resistance"] for row in uniformity_rows}):
        draw_uniformity(uniformity_rows, resistance, args.workdir)


if __name__ == "__main__":
    main()
