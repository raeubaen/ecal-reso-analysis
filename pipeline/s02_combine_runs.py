#!/usr/bin/env python3
"""
Stage 2 -- combine the runs of each energy, and save.

Reads 01_dcb_per_run.csv. For every (resistance, energy, variation) the per-run values
of sigma/mu [%] are combined into one point. The weights follow the recipe of the
selection, which is what the two flat chains do:

  codiceA  (hodoscope)  weights 1/err^2 (resolution_hodo.py)
  uniforme (centroid)   weights = number of selected events of the run (uniformita_pos.py)

and for both:
  stat         max(error of the weighted mean, weighted scatter of the values): the
               scatter already contains the fit noise and the run-to-run spread, and is
               undefined with one run
  drift        the extra error that brings the per-run sigma/mu to chi2/ndf = 1 against
               a constant (drift_dcb_all.syst_for_unit_chi2 on sigma/mu, as both flat
               scripts do after merge #2); 0 with one run
  syst_tails   |mean with free tails - mean with the tails held at the pooled values|.
               resolution_hodo.py computes it and does not write it anywhere: it is
               saved here for information and does NOT enter any error bar
  vtx_syst     codiceA only, on the nominal row: sqrt(sum over the four shifted windows
               of (mean_variation - mean_nominal)^2 / 4)
  pooled       codiceA only: where no run reaches 300 events the single pooled fit is
               the point (run = 0, pooled = 1)

Writes:
  02_per_energy.csv   one row per (resistance, energy, variation)
  02_per_run.csv      the nominal variation run by run, with the deviation and the
                      pull of peak and sigma from the energy mean (which run pulls the
                      drift)
  runs/*.png          peak, sigma and sigma/mu against the run, one canvas per energy
  drift_check_<R>ohm.png   per-run sigma/mu against the constant, with the drift band
"""

import argparse
import math
import os
from collections import defaultdict

import numpy as np
import ROOT

import common

ENERGY_COLUMNS = ("resistance", "energy", "energy_true", "selection", "recipe", "variation",
                  "window", "n_run", "pooled", "n_events", "sigma_over_mu", "err_mean",
                  "scatter", "stat", "syst_tails", "drift", "chi2_ndf_drift", "vtx_syst",
                  "peak_mean", "runs")
RUN_COLUMNS = ("resistance", "energy", "energy_true", "selection", "run", "n_selected",
               "n_events", "weight", "sigma_pct", "err_sigma_pct", "sigma_fixed_tails",
               "d_tails", "peak", "err_peak", "peak_dev_pct", "peak_pull", "sigma_dev_pct",
               "sigma_pull", "E_sigma_mean", "E_drift_pct", "E_chi2_drift", "E_n_run")


def weights_for(recipe, rows):
    if recipe == "codiceA":
        return [1. / row["err_sigma_over_mu"] ** 2 if row["err_sigma_over_mu"] > 0 else 0.
                for row in rows]
    return [float(row["n_selected"]) for row in rows]


def combine(point_rows, recipe, nominal_tails):
    """One (resistance, energy, variation): returns the summary dict and the run rows."""
    runs = sorted(row["run"] for row in point_rows
                  if row["run"] > 0 and row["tails_mode"] == nominal_tails)
    by_run = {row["run"]: row for row in point_rows
              if row["run"] > 0 and row["tails_mode"] == nominal_tails}
    fixed_by_run = {row["run"]: row for row in point_rows
                    if row["run"] > 0 and row["tails_mode"] == "fixed"}
    pooled_rows = [row for row in point_rows if row["run"] == 0]
    nominal_rows = [by_run[run] for run in runs]
    pooled = 0
    if not nominal_rows and recipe == "codiceA" and pooled_rows:
        nominal_rows, pooled = pooled_rows[:1], 1
    if not nominal_rows:
        return None, []

    values = [row["sigma_over_mu"] for row in nominal_rows]
    errors = [row["err_sigma_over_mu"] for row in nominal_rows]
    weights = weights_for(recipe, nominal_rows)
    mean, err_mean = common.weighted_mean(values, errors, weights)
    scatter = common.weighted_scatter(values, weights)
    stat = max(err_mean, scatter) if np.isfinite(scatter) else err_mean

    syst_tails = 0.
    fixed_runs = [run for run in runs if run in fixed_by_run]
    if fixed_runs and nominal_tails == "free":
        fixed_values = [fixed_by_run[run]["sigma_over_mu"] for run in fixed_runs]
        free_errors = [by_run[run]["err_sigma_over_mu"] for run in fixed_runs]
        fixed_weights = weights_for(recipe, [by_run[run] for run in fixed_runs])
        mean_fixed = common.weighted_mean(fixed_values, free_errors, fixed_weights)[0]
        if np.isfinite(mean_fixed):
            syst_tails = abs(mean - mean_fixed)

    drift, chi2_ndf = 0., np.nan
    if len(values) > 1:
        drift, chi2_ndf, _mean_with_drift = common.syst_for_unit_chi2(values, errors)

    peaks = np.array([row["peak"] for row in nominal_rows])
    err_peaks = np.array([row["err_peak"] for row in nominal_rows])
    peak_weights = 1. / np.maximum(err_peaks, 1e-12) ** 2
    peak_mean = float((peaks * peak_weights).sum() / peak_weights.sum())

    summary = dict(n_run=0 if pooled else len(nominal_rows), pooled=pooled,
                   n_events=int(sum(row["n_selected"] for row in nominal_rows)),
                   sigma_over_mu=mean, err_mean=err_mean, scatter=scatter, stat=stat,
                   syst_tails=syst_tails, drift=drift, chi2_ndf_drift=chi2_ndf, vtx_syst=0.,
                   peak_mean=peak_mean,
                   runs=";".join(str(row["run"]) for row in nominal_rows))
    run_rows = []
    for row, weight in zip(nominal_rows, weights):
        fixed = fixed_by_run.get(row["run"])
        sigma_fixed = fixed["sigma_over_mu"] if fixed else np.nan
        run_rows.append(dict(
            run=row["run"], n_selected=row["n_selected"], n_events=row["n_events"],
            weight=weight, sigma_pct=row["sigma_over_mu"], err_sigma_pct=row["err_sigma_over_mu"],
            sigma_fixed_tails=sigma_fixed, d_tails=abs(row["sigma_over_mu"] - sigma_fixed),
            peak=row["peak"], err_peak=row["err_peak"],
            peak_dev_pct=100. * (row["peak"] / peak_mean - 1.),
            peak_pull=(row["peak"] - peak_mean) / row["err_peak"] if row["err_peak"] > 0 else np.nan,
            sigma_dev_pct=100. * (row["sigma_over_mu"] / mean - 1.) if mean > 0 else np.nan,
            sigma_pull=((row["sigma_over_mu"] - mean) / row["err_sigma_over_mu"]
                        if row["err_sigma_over_mu"] > 0 else np.nan),
            E_sigma_mean=mean, E_drift_pct=drift, E_chi2_drift=chi2_ndf,
            E_n_run=summary["n_run"]))
    return summary, run_rows


def draw_runs(resistance, energy, run_rows, summary, pooled_row, outdir, selection):
    """peak, sigma and sigma/mu against the run, with the pooled fit as reference."""
    n_runs = len(run_rows)
    x_values = np.arange(n_runs, dtype=float)
    labels = [row["run"] for row in run_rows]
    canvas = ROOT.TCanvas(f"runs_{resistance}_{energy}", "", max(800, 120 * n_runs), 900)
    canvas.Divide(1, 3)
    panels = (("peak [ADC]", [row["peak"] for row in run_rows], [row["err_peak"] for row in run_rows],
               pooled_row["peak"] if pooled_row else None),
              ("#sigma/#mu [%]", [row["sigma_pct"] for row in run_rows],
               [row["err_sigma_pct"] for row in run_rows],
               pooled_row["sigma_over_mu"] if pooled_row else None),
              ("#sigma/#mu / mean - 1 [%]", [row["sigma_dev_pct"] for row in run_rows],
               [100 * row["err_sigma_pct"] / summary["sigma_over_mu"] for row in run_rows], 0.))
    for index, (label, y_values, y_errors, reference) in enumerate(panels):
        pad = canvas.cd(index + 1)
        pad.SetLeftMargin(0.12)
        pad.SetBottomMargin(0.22)
        values, errors = np.array(y_values), np.array(y_errors)
        low = min(np.min(values - errors), reference if reference is not None else np.inf)
        high = max(np.max(values + errors), reference if reference is not None else -np.inf)
        margin = 0.15 * (high - low) + 1e-9
        frame = common.run_axis_frame(labels, low - margin, high + margin,
                                      f"{resistance} #Omega, {energy} GeV, {selection};run;{label}")
        frame.Draw()
        if index == 2 and summary["drift"] > 0:
            drift_rel = 100 * summary["drift"] / summary["sigma_over_mu"]
            band = common.keep(ROOT.TBox(-0.5, -drift_rel, n_runs - 0.5, drift_rel))
            band.SetFillColorAlpha(ROOT.kRed, 0.15)
            band.Draw()
        if reference is not None:
            line = common.keep(ROOT.TLine(-0.5, reference, n_runs - 0.5, reference))
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray + 2)
            line.Draw()
        graph = common.keep(common.make_graph(x_values, values, errors))
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(common.COLOUR[resistance])
        graph.SetLineColor(common.COLOUR[resistance])
        graph.Draw("P")
        if index == 2:
            chi2 = summary["chi2_ndf_drift"]
            note = (f"#chi^{{2}}/ndf {chi2:.2f} ({n_runs - 1} ndf)  #rightarrow  drift "
                    f"{summary['drift']:.4f} % of #sigma/#mu" if np.isfinite(chi2)
                    else "one run: drift undefined")
            text = common.keep(ROOT.TLatex(0.14, 0.85, note))
            text.SetNDC()
            text.SetTextSize(0.05)
            text.Draw()
    common.save_canvas(canvas, os.path.join(outdir, "runs",
                                            f"peak_sigma_vs_run_{energy}GeV_{resistance}ohm.png"))


def draw_drift_check(resistance, energy_points, outdir, selection):
    """One pad per energy with more than one run: sigma/mu of each run relative to the
    weighted mean, the constant and the drift band that brings chi2/ndf to 1."""
    canvas, _rows = common.canvas_grid(f"drift_{resistance}", len(energy_points), columns=4)
    for index, (energy, summary, run_rows) in enumerate(energy_points):
        pad = canvas.cd(index + 1)
        pad.SetBottomMargin(0.2)
        n_runs = len(run_rows)
        values = np.array([row["sigma_dev_pct"] for row in run_rows])
        errors = np.array([100 * row["err_sigma_pct"] / summary["sigma_over_mu"] for row in run_rows])
        drift_rel = 100 * summary["drift"] / summary["sigma_over_mu"]
        chi2 = summary["chi2_ndf_drift"]
        verdict = (f"#chi^{{2}}/ndf {chi2:.2f} > 1 #rightarrow drift {summary['drift']:.4f} %"
                   if summary["drift"] > 0 else f"#chi^{{2}}/ndf {chi2:.2f} #leq 1 #rightarrow drift 0")
        low = min(np.min(values - errors), -1.2 * drift_rel)
        high = max(np.max(values + errors), 1.2 * drift_rel)
        margin = 0.1 * (high - low) + 1e-9
        frame = common.run_axis_frame([row["run"] for row in run_rows], low - margin, high + margin,
                                      f"{energy} GeV: {verdict};run;#sigma/#mu / mean - 1  [%]")
        frame.Draw()
        if summary["drift"] > 0:
            band = common.keep(ROOT.TBox(-0.5, -drift_rel, n_runs - 0.5, drift_rel))
            band.SetFillColorAlpha(ROOT.kRed, 0.15)
            band.Draw()
        line = common.keep(ROOT.TLine(-0.5, 0., n_runs - 0.5, 0.))
        line.SetLineColor(ROOT.kRed)
        line.Draw()
        graph = common.keep(common.make_graph(np.arange(n_runs, dtype=float), values, errors))
        graph.SetMarkerStyle(20)
        graph.SetMarkerColor(ROOT.kAzure + 2)
        graph.SetLineColor(ROOT.kAzure + 2)
        graph.Draw("P")
    canvas.cd()
    header = common.keep(ROOT.TLatex(0.01, 0.985, f"{resistance} #Omega, {selection}: per-run #sigma/#mu against a constant"))
    header.SetNDC()
    header.SetTextSize(0.018)
    header.Draw()
    common.save_canvas(canvas, os.path.join(outdir, f"drift_check_{resistance}ohm.png"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True, help="the --outdir of stage 1")
    parser.add_argument("--nominal-tails", choices=("free", "fixed"), default="free",
                        help="which per-run fit is the nominal one (free, as --tails both)")
    args = parser.parse_args()
    common.style()

    fit_rows = common.read_csv(common.require(os.path.join(args.workdir, "01_dcb_per_run.csv")))
    windows = {(row["resistance"], row["energy"]): row
               for row in common.read_csv(common.require(os.path.join(args.workdir, "01_windows.csv")))}
    grouped = defaultdict(list)
    for row in fit_rows:
        grouped[(row["resistance"], row["energy"], row["variation"])].append(row)

    energy_rows, run_rows_all = [], []
    summaries = {}
    for (resistance, energy, variation), point_rows in sorted(grouped.items()):
        recipe, selection = point_rows[0]["recipe"], point_rows[0]["selection"]
        summary, run_rows = combine(point_rows, recipe, args.nominal_tails)
        if summary is None:
            print(f"  {resistance} ohm {energy:>4} GeV {variation}: no fit, point dropped")
            continue
        summary.update(resistance=resistance, energy=energy, energy_true=common.true_energy(energy),
                       selection=selection, recipe=recipe, variation=variation,
                       window=windows[(resistance, energy)]["window"])
        if variation == "nominal":
            summary["n_events"] = windows[(resistance, energy)]["n_selected"]
        summaries[(resistance, energy, variation)] = summary
        energy_rows.append(summary)
        if variation == "nominal":
            for run_row in run_rows:
                run_row.update(resistance=resistance, energy=energy,
                               energy_true=common.true_energy(energy), selection=selection)
            run_rows_all += run_rows
            pooled = [row for row in point_rows if row["run"] == 0 and row["tails_mode"] == "free"]
            if run_rows and not summary["pooled"]:
                draw_runs(resistance, energy, run_rows, summary, pooled[0] if pooled else None,
                          args.workdir, selection)

    # vertex systematic of the hodoscope window, on the nominal row
    for (resistance, energy, variation), summary in summaries.items():
        if variation != "nominal" or summary["recipe"] != "codiceA":
            continue
        shifted = [summaries[key]["sigma_over_mu"] for key in summaries
                   if key[:2] == (resistance, energy) and key[2] != "nominal"]
        deviations = [(value - summary["sigma_over_mu"]) ** 2 for value in shifted if np.isfinite(value)]
        summary["vtx_syst"] = math.sqrt(sum(deviations) / 4.) if deviations else 0.
        print(f"  {resistance} ohm {energy:>4} GeV: sigma/mu {summary['sigma_over_mu']:.4f} % "
              f"(stat {summary['stat']:.4f}, drift {summary['drift']:.4f}, "
              f"vtx {summary['vtx_syst']:.4f}, {summary['n_run']} run)")

    energy_rows.sort(key=lambda row: (row["resistance"], row["energy"], row["variation"] != "nominal",
                                      row["variation"]))
    common.write_csv(os.path.join(args.workdir, "02_per_energy.csv"), energy_rows, ENERGY_COLUMNS)
    common.write_csv(os.path.join(args.workdir, "02_per_run.csv"), run_rows_all, RUN_COLUMNS)

    for resistance in sorted({row["resistance"] for row in energy_rows}):
        points = []
        for (this_resistance, energy, variation), summary in sorted(summaries.items()):
            if this_resistance != resistance or variation != "nominal" or summary["n_run"] < 2:
                continue
            runs = [row for row in run_rows_all
                    if row["resistance"] == resistance and row["energy"] == energy]
            points.append((energy, summary, runs))
        if points:
            draw_drift_check(resistance, points, args.workdir, points[0][1]["selection"])


if __name__ == "__main__":
    main()
