#!/usr/bin/env python3
"""
Stage 1 -- double Crystal Ball fit of the amplitude, run by run, and save.

  hodoscope   the window on the hodoscope of resolution_hodo.py: vertex of the response
              parabola +- 4mm in x and in y, plus the four windows shifted by
              +- 1 mm that feed the vertex systematic (variations x_low, x_high, y_low,
              y_high). 

For every (resistance, energy, variation):
  * one POOLED fit of all the runs together (run = 0), whose tails seed the "fixed"
    model;
  * one fit per run with the tails free and, with --tails both, one with the tails
    held at the pooled values.

Writes, in --outdir:
  01_dcb_per_run.csv   one row per fit: peak, sigma, their HESSE errors, sigma/mu in %,
                       chi2, ndf, tails, window, number of events
  01_windows.csv       one row per (resistance, energy): the position window used and,
                       for the hodoscope, the parabola-scan diagnostics
  01_dcb_fits.root     every histogram with its fitted TF1
  dcb/*.png            one panel per (resistance, energy): the per-run fits + pooled

Usage:
  python3 s01_fit_dcb_per_run.py --base <dir> --outdir out/hodoscope\\
"""

import argparse
import os

import numpy as np
import ROOT

import common
import hodoscope_window
from common import runsets

FIT_COLUMNS = ("resistance", "energy", "energy_true", "amplitude",
               "variation", "run", "tails_mode", "n_selected", "n_events", "n_bins",
               "peak", "err_peak", "sigma", "err_sigma", "sigma_over_mu", "err_sigma_over_mu",
               "chi2", "ndf", "alpha_l", "alpha_h", "n_l", "n_h", "window_lo", "window_hi",
               "valid")
WINDOW_COLUMNS = ("resistance", "energy", "energy_true", "window",
                  "n_base", "n_selected", "skipped", "reason",
                  "x_lo", "x_hi", "y_lo", "y_hi", "x_vertex", "y_vertex", "x_width", "y_width",
                  "x_ok", "y_ok", "x_why", "y_why", "fallback")


def fit_row(fit, resistance, energy, args, variation, run, tails_mode, n_selected):
    value, error = common.relative_width(fit)
    row = dict(resistance=resistance, energy=energy, energy_true=common.true_energy(energy),
               amplitude=args.amplitude, variation=variation, run=run, tails_mode=tails_mode,
               n_selected=int(n_selected), n_events=fit["n_events"], n_bins=fit["n_bins"],
               peak=fit["peak"], err_peak=fit["err_peak"], sigma=fit["sigma"],
               err_sigma=fit["err_sigma"], sigma_over_mu=value, err_sigma_over_mu=error,
               chi2=fit["chi2"], ndf=fit["ndf"], window_lo=fit["lo"], window_hi=fit["hi"],
               valid=int(fit["valid"]))
    row.update(fit["tails"])
    return row


def store_fit(root_file, fit, resistance, energy, variation, run, tails_mode):
    name = f"{resistance}_{energy}_{variation}_{run}_{tails_mode}"
    root_file.cd()
    fit["histogram"].SetName("h_" + name)
    fit["histogram"].SetTitle(f"{resistance} ohm {energy} GeV {variation} run {run} {tails_mode}")
    fit["histogram"].Write()
    fit["function"].SetName("f_" + name)
    fit["function"].Write()


def draw_fit_panel(fits_by_run, pooled_fit, resistance, energy, args):
    """The per-run fits of the nominal selection on one canvas, the pooled fit last."""
    labels = [str(run) for run in fits_by_run] + ["ALL RUNS"]
    fits = list(fits_by_run.values()) + [pooled_fit]
    canvas, _rows = common.canvas_grid(f"dcb_{resistance}_{energy}", len(labels))
    for index, (label, fit) in enumerate(zip(labels, fits)):
        pad = canvas.cd(index + 1)
        pad.SetLeftMargin(0.13)
        if fit is None:
            text = common.keep(ROOT.TLatex(0.3, 0.5, f"{label}: fit failed"))
            text.SetNDC()
            text.Draw()
            continue
        histogram, function = fit["histogram"], fit["function"]
        histogram.SetMarkerStyle(20)
        histogram.SetMarkerSize(0.5)
        histogram.SetLineColor(ROOT.kBlack)
        histogram.GetXaxis().SetTitle("amplitude [ADC]")
        histogram.SetTitle(f"{label}: #mu = {fit['peak']:.1f}  #sigma = {fit['sigma']:.2f}"
                           f"  #chi^{{2}}/ndf = {fit['chi2'] / fit['ndf']:.2f}")
        histogram.Draw("E1")
        function.SetLineColor(ROOT.kRed)
        function.SetLineWidth(2)
        function.SetNpx(600)
        function.Draw("same")
    title = f"{resistance} #Omega, {energy} GeV"
    canvas.cd()
    header = common.keep(ROOT.TLatex(0.01, 0.985, title))
    header.SetNDC()
    header.SetTextSize(0.018)
    header.Draw()
    common.save_canvas(canvas, os.path.join(args.outdir, "dcb",
                                            f"dcb_fits_{energy}GeV_{resistance}ohm.png"))
    common.save_canvas(canvas, os.path.join(args.outdir, "dcb",
                                            f"dcb_fits_{energy}GeV_{resistance}ohm.pdf"))
    common.save_canvas(canvas, os.path.join(args.outdir, "dcb",
                                            f"dcb_fits_{energy}GeV_{resistance}ohm.root"))



def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", required=True, help="directory containing reco_<R>ohm/")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--amplitude", choices=("a3x3", "atot", "a1x1", "a5x5"), default="a3x3",
                        help="a3x3 = sum of the 3x3 matrix rebuilt from A (default, as "
                             "resolution_hodo.py); atot = the A_tot branch")
    parser.add_argument("--resistances", nargs="+", type=int, default=[340, 400, 500])
    parser.add_argument("--eta-center", type=int, default=18)
    parser.add_argument("--phi-center", type=int, default=6)
    parser.add_argument("--yplane", choices=("y1", "y2"), default="y1")
    parser.add_argument("--tails", choices=("free", "fixed", "both"), default="both",
                        help="DCB tails per run: free, held at the pooled values, or both")

    parser.add_argument("--half", required=False, type=float, default=4, help="hodoscope selection half window, in mm")
    parser.add_argument("--exclude-runs", nargs="*", type=int, default=[])
    parser.add_argument("--exclude", nargs="*", default=["340:275"],
                        help="R:E points dropped entirely (default 340:275, as both drivers)")
    parser.add_argument("--fallback-file", type=str, default="fallback_hodo_26.py")

    runsets.add_argument(parser)
    args = parser.parse_args()

    half = args.half
    dropped, kept_only = runsets.resolve(args.runset, args.exclude_runs)
    print(dropped, kept_only)
    excluded_points = common.parse_excluded_points(args.exclude)
    os.makedirs(args.outdir, exist_ok=True)
    common.style()
    root_file = ROOT.TFile(os.path.join(args.outdir, "01_dcb_fits.root"), "RECREATE")

    fit_rows, window_rows = [], []
    print("pre for")
    for resistance, energy, path in common.resistance_energy_pairs(args.base, args.resistances, excluded_points):

        print(f"[{resistance} ohm {energy:>4} GeV] {os.path.basename(path)}", flush=True)
        events = common.read_events(path, args.eta_center, args.phi_center, args.amplitude)
        base = (events["A_tot"] > common.A_TOT_MIN) & common.runset_mask(events["run"], dropped,
                                                                          kept_only)

        window_row = dict(resistance=resistance, energy=energy,
                          energy_true=common.true_energy(energy),
                          n_base=int(base.sum()), skipped=0, reason="",
                          fallback="")

        hodo_x, hodo_y = common.hodoscope_xy(events, args.yplane)

        info = hodoscope_window.hodoscope_windows(hodo_x, hodo_y, events["A_tot"], base, resistance, energy, half, args.outdir, args.fallback_file)

        for coordinate in ("x", "y"):
            scan = info["scan"][coordinate]
            window_row.update({f"{coordinate}_vertex": scan["vertex"],
                               f"{coordinate}_width": scan["width"],
                               f"{coordinate}_ok": int(scan["ok"]),
                               f"{coordinate}_why": info["why"][coordinate]})
        window_row["fallback"] = "+".join(info["fallback"])
        window_row["window"] = hodoscope_window.window_label(info["fallback"])
        for coordinate in info["fallback"]:
            reason = info["why"][coordinate]
            print(f"    {coordinate}: hand-set vertex of resolution_hodo.py"
                  + (f" (scan failed: {reason})" if reason else " (overrides a successful scan)"))
        missing = [c for c in ("x", "y") if info["windows"][c] is None]
        if missing:
            reason = "; ".join(f"{c}: {info['why'][c]}" for c in missing)
            print(f"    SKIPPED, no window in {'+'.join(missing)} ({reason})")
            window_row.update(skipped=1, reason=reason, n_selected=0)
            window_rows.append(window_row)
            continue
        window_x, window_y = info["windows"]["x"], info["windows"]["y"]
        window_row.update(x_lo=window_x[0], x_hi=window_x[1], y_lo=window_y[0],
                          y_hi=window_y[1])
        cuts = hodoscope_window.window_masks(hodo_x, hodo_y, base, window_x, window_y)
        if cuts["nominal"].sum() < common.MIN_EVENTS_POOLED:
            print(f"    SKIPPED, only {cuts['nominal'].sum()} events after the hodoscope cut")
            window_row.update(skipped=1, n_selected=int(cuts["nominal"].sum()),
                              reason="fewer than 500 events in the window")
            window_rows.append(window_row)
            continue

        window_row["n_selected"] = int(cuts["nominal"].sum())
        window_rows.append(window_row)

        amplitude, run = events["amplitude"], events["run"]
        nominal_fits, pooled_nominal = {}, None
        for variation, cut in cuts.items():
            pooled_fit = None
            if cut.sum() >= common.MIN_EVENTS_POOLED:
                pooled_fit = common.fit_dcb(amplitude[cut], energy, resistance, args.amplitude)
            tails = pooled_fit["tails"] if (pooled_fit and args.tails != "free") else None
            if pooled_fit is not None:
                fit_rows.append(fit_row(pooled_fit, resistance, energy, args, variation, 0,
                                        "free", cut.sum()))
                store_fit(root_file, pooled_fit, resistance, energy, variation, 0, "free")
            for this_run in sorted(int(value) for value in np.unique(run[cut])):
                in_run = cut & (run == this_run)
                if in_run.sum() < common.MIN_EVENTS_PER_RUN:
                    print("in_run.sum() < common.MIN_EVENTS_PER_RUN")
                    continue
                values = amplitude[in_run]
                print(values)
                if args.tails == "fixed":
                    modes = [("fixed", tails)] if tails else []
                elif args.tails == "both" and tails:
                    modes = [("free", None), ("fixed", tails)]
                else:
                    modes = [("free", None)]
                for tails_mode, fixed_tails in modes:
                    fit = common.fit_dcb(values, energy, resistance, args.amplitude, fixed_tails=fixed_tails)
                    if fit is None:
                        print("BAD FIT!")
                        continue
                    fit_rows.append(fit_row(fit, resistance, energy, args, variation, this_run,
                                            tails_mode, in_run.sum()))
                    store_fit(root_file, fit, resistance, energy, variation, this_run, tails_mode)
                    if variation == "nominal" and tails_mode == modes[0][0]:
                        nominal_fits[this_run] = fit
            if variation == "nominal":
                pooled_nominal = pooled_fit
        n_runs = len(nominal_fits)
        print(f"    {n_runs} run fitted, {cuts['nominal'].sum()} events in the nominal cut"
              + (f", pooled sigma/mu {common.relative_width(pooled_nominal)[0]:.4f} %"
                 if pooled_nominal else ""), flush=True)
        if nominal_fits or pooled_nominal:
            draw_fit_panel(nominal_fits, pooled_nominal, resistance, energy, args)

    root_file.Close()
    common.write_csv(os.path.join(args.outdir, "01_dcb_per_run.csv"), fit_rows, FIT_COLUMNS)
    common.write_csv(os.path.join(args.outdir, "01_windows.csv"), window_rows, WINDOW_COLUMNS)


if __name__ == "__main__":
    main()
