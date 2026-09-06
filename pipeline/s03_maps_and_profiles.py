#!/usr/bin/env python3
"""
Stage 3 -- 2D maps and response profiles against the centroid, and save.

Centroid pass only (the hodoscope chain does not correct for the position). Reads the
events and the per-run fits of stage 1 (peak, sigma and fit window of every run) and
produces three tables and the plots.

  03_moments.csv     per run: the sufficient statistics of the response surface
                     f = a0 + a1 u + a2 v + a3 u^2 + a4 v^2 + a5 uv on the events inside
                     the fit window of the run, normalised to the peak of the run:
                     M = sum x x' (6x6) and b = sum x a, with x = (1, u, v, u^2, v^2, uv).
                     Stages 4 and 5 solve and average the surfaces from these numbers
                     without reading an event again (uniformita_maps.py, stage moments).
  03_profiles.csv    <A> against pos_eta - 18 (|pos_phi - 6| <= 0.2) and against
                     pos_phi - 6, every event normalised to the peak of its own run, with
                     the quadratic a + b x + c x^2 over |x| <= 0.3 and its relative
                     curvature c/a in %/crystal^2 (profili_pernorm.py). Same columns as
                     profili_pernorm.csv: this is the curvature the hodoscope window needs.
  03_parabola_centroide.csv   the same profiles on the pooled events with the parabola
                     p1 + p2 (x - p0)^2 of drift_dcb_all.py (vertex p0, maximum p1,
                     curvature p2), fitted with MIGRAD

  maps/*.png         occupancy and <A> map in (pos_eta, pos_phi), the centroid of each
                     run, the two profiles with the parabola (drift_dcb_all.py figures)
  profiles/*.png     the per-run-normalised profiles with the quadratic, the curvature
                     against the energy, all the parabolas of a resistance together

The quadratic fits are weighted linear least squares (ROOT linear fitter, "pol2").
"""

import argparse
import os
from collections import defaultdict

import numpy as np
import ROOT

import common
from common import runsets

PROFILE_HALF, PROFILE_NBINS, FIT_HALF, PROFILE_MIN_PER_BIN = 0.6, 96, 0.3, 20
N_SIGMA_WINDOW = 10.                 # |A - peak_run| < 10 sigma_run for the profiles
MIN_EVENTS_MOMENTS = 100             # a run needs this many events in the window for M, b
MIN_EVENTS_RUN_CURVE = 200           # ... and this many to be drawn as its own curve
MAP_HALF, MAP_NBINS, MAP_MIN_PER_BIN = 0.6, 200, 3
PALETTE_HALF = 0.03                  # colour scale = inclusive mean +- 3 %
# noise threshold of drift_dcb_all.py for the maps: 80 ADC up to 50 GeV, 200 above,
# or 5 % of the nominal amplitude when larger
A_BASE_LOW, A_BASE_HIGH, E_LOW, A_FRACTION = 80., 200., 50., 0.05
COORDINATES = (("pos_eta", "u", "v", "pos_eta - 18", "|pos_phi - 6| <= 0.2"),
               ("pos_phi", "v", "u", "pos_phi - 6", "|pos_eta - 18| <= 0.2"))

MOMENT_COLUMNS = (["resistance", "energy", "run", "n_selected", "n_window", "peak", "sigma",
                   "sigma_over_mu", "err_sigma_over_mu", "window_lo", "window_hi"]
                  + [f"M_{row}{column}" for row in range(6) for column in range(6)]
                  + [f"b_{index}" for index in range(6)])
PROFILE_COLUMNS = ("resistance,energy,coord,nrun,nev,a,b,c,err_a,err_b,err_c,"
                   "rel_pct,err_rel_pct,chi2,ndf,rel_pct_raw,chi2_raw,ndf_raw").split(",")
PARABOLA_COLUMNS = ("resistance", "energy", "coord", "n_bin", "p0_vertice", "err_p0",
                    "p1_massimo", "err_p1", "p2_curvatura", "err_p2", "chi2", "ndf")


def map_threshold(energy, resistance):
    base = A_BASE_LOW if energy <= E_LOW else A_BASE_HIGH
    fraction = A_FRACTION * common.SCALE[resistance] * energy
    return max(base, fraction) if fraction > A_BASE_LOW else base


def design(u_values, v_values):
    return np.column_stack([np.ones_like(u_values), u_values, v_values,
                            u_values * u_values, v_values * v_values, u_values * v_values])


def binned_profile(coordinate, response):
    """Mean and error on the mean of the response in PROFILE_NBINS bins over
    +- PROFILE_HALF; only bins with at least PROFILE_MIN_PER_BIN events."""
    counts, edges = np.histogram(coordinate, bins=PROFILE_NBINS, range=(-PROFILE_HALF, PROFILE_HALF))
    sums, _ = np.histogram(coordinate, bins=PROFILE_NBINS, range=(-PROFILE_HALF, PROFILE_HALF),
                           weights=response)
    squares, _ = np.histogram(coordinate, bins=PROFILE_NBINS, range=(-PROFILE_HALF, PROFILE_HALF),
                              weights=response ** 2)
    centres = 0.5 * (edges[:-1] + edges[1:])
    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / np.maximum(counts, 1)
        variances = squares / np.maximum(counts, 1) - means ** 2
        errors = np.sqrt(np.maximum(variances, 0) / np.maximum(counts, 1))
    good = (counts >= PROFILE_MIN_PER_BIN) & (errors > 0)
    return centres[good], means[good], errors[good], int(counts[good].sum())


def quadratic_fit(centres, means, errors, half=FIT_HALF):
    """a + b x + c x^2 over |x| <= half, weighted linear least squares."""
    inside = np.abs(centres) <= half
    if inside.sum() < 8:
        return None
    graph = common.make_graph(centres[inside], means[inside], errors[inside])
    quadratic = ROOT.TF1(common.unique_name("prof_pol2"), "pol2", -half, half)
    graph.Fit(quadratic, "QN")
    coefficients = [quadratic.GetParameter(index) for index in range(3)]
    coefficient_errors = [quadratic.GetParError(index) for index in range(3)]
    constant, curvature = coefficients[0], coefficients[2]
    return dict(coefficients=coefficients, errors=coefficient_errors,
                chi2=quadratic.GetChisquare(), ndf=int(inside.sum()) - 3,
                rel=100 * curvature / constant, err_rel=100 * coefficient_errors[2] / abs(constant),
                function=quadratic, graph=graph)


def per_run_fits(fit_rows, resistance, energy):
    """{run: row} of the nominal free per-run fits of stage 1."""
    return {row["run"]: row for row in fit_rows
            if row["resistance"] == resistance and row["energy"] == energy
            and row["variation"] == "nominal" and row["tails_mode"] == "free" and row["run"] > 0}


# ------------------------------------------------------------------ moments
def moments_rows(events, selected, fits, resistance, energy):
    rows = []
    run = events["run"]
    for this_run, fit in sorted(fits.items()):
        in_run = selected & (run == this_run)
        in_window = in_run & (events["amplitude"] >= fit["window_lo"]) & (events["amplitude"] <= fit["window_hi"])
        if in_window.sum() < MIN_EVENTS_MOMENTS:
            print(f"    run {this_run}: only {in_window.sum()} events in the fit window, no surface")
            continue
        matrix = design(events["u"][in_window], events["v"][in_window])
        normalised = events["amplitude"][in_window] / fit["peak"]
        moment_matrix = matrix.T @ matrix
        moment_vector = matrix.T @ normalised
        row = dict(resistance=resistance, energy=energy, run=this_run, n_selected=int(in_run.sum()),
                   n_window=int(in_window.sum()), peak=fit["peak"], sigma=fit["sigma"],
                   sigma_over_mu=fit["sigma_over_mu"], err_sigma_over_mu=fit["err_sigma_over_mu"],
                   window_lo=fit["window_lo"], window_hi=fit["window_hi"])
        for row_index in range(6):
            for column_index in range(6):
                row[f"M_{row_index}{column_index}"] = float(moment_matrix[row_index, column_index])
            row[f"b_{row_index}"] = float(moment_vector[row_index])
        rows.append(row)
    return rows


# ------------------------------------------------------------------ normalised profiles
def normalised_profiles(events, base, fits, resistance, energy, outdir):
    """profili_pernorm.analyse: every event scaled to the peak of its run, profile in
    each coordinate, quadratic fit, one canvas."""
    run, amplitude = events["run"], events["amplitude"]
    runs = sorted(this_run for this_run in fits if (base & (run == this_run)).any())
    if not runs:
        return []
    peak_of_event = np.ones(len(amplitude))
    sigma_of_event = np.ones(len(amplitude))
    for this_run in runs:
        in_run = run == this_run
        peak_of_event[in_run] = fits[this_run]["peak"]
        sigma_of_event[in_run] = fits[this_run]["sigma"]
    keep = base & np.isin(run, runs)
    in_window = keep & (np.abs(amplitude - peak_of_event) < N_SIGMA_WINDOW * sigma_of_event)
    events_per_run = np.array([(in_window & (run == this_run)).sum() for this_run in runs], float)
    reference_peak = float((events_per_run * [fits[this_run]["peak"] for this_run in runs]).sum()
                           / events_per_run.sum())
    normalised = amplitude / peak_of_event * reference_peak

    canvas = ROOT.TCanvas(f"profiles_{resistance}_{energy}", "", 1500, 560)
    canvas.Divide(2, 1)
    rows = []
    for pad_index, (name, own, other, label, other_cut) in enumerate(COORDINATES):
        pad = canvas.cd(pad_index + 1)
        pad.SetLeftMargin(0.12)
        pad.SetBottomMargin(0.12)
        coordinate, other_coordinate = events[own], events[other]
        in_profile = in_window & (np.abs(other_coordinate) <= common.HALF_WINDOW["centroid"])
        centres, means, errors, n_events = binned_profile(coordinate[in_profile], normalised[in_profile])
        fit = quadratic_fit(centres, means, errors)
        raw_fit = quadratic_fit(*binned_profile(coordinate[in_profile], amplitude[in_profile])[:3])
        frame = common.keep(ROOT.TH2F(common.unique_name("frame"), "", 10, -PROFILE_HALF, PROFILE_HALF,
                                      10, means.min() * 0.97 if len(means) else 0,
                                      means.max() * 1.03 if len(means) else 1))
        frame.SetTitle(f"{resistance} #Omega, {energy} GeV: {other_cut}, |A - #mu_{{run}}| < "
                       f"{N_SIGMA_WINDOW:.0f}#sigma_{{run}}, {len(runs)} runs, N = {int(in_window.sum())}"
                       f";{label} [crystal units];#LTA A / peak_{{run}} #GT #times ref [ADC]")
        frame.Draw()
        legend = common.keep(ROOT.TLegend(0.62, 0.13, 0.89, 0.13 + 0.035 * min(len(runs) + 1, 9)))
        legend.SetTextSize(0.028)
        legend.SetBorderSize(0)
        for run_index, this_run in enumerate(runs):
            in_run_profile = in_profile & (run == this_run)
            if in_run_profile.sum() < MIN_EVENTS_RUN_CURVE:
                continue
            run_centres, run_means, _errors, _n = binned_profile(coordinate[in_run_profile],
                                                                  normalised[in_run_profile])
            if len(run_centres) < 2:
                continue
            curve = common.keep(common.make_graph(run_centres, run_means))
            curve.SetLineColor(ROOT.TColor.GetPalette()[int(240 * run_index / max(len(runs) - 1, 1))])
            curve.SetLineWidth(1)
            curve.Draw("L")
            if len(runs) <= 8:
                legend.AddEntry(curve, f"run {this_run}", "l")
        if len(centres):
            points = common.keep(common.make_graph(centres, means, errors))
            points.SetMarkerStyle(20)
            points.SetMarkerSize(0.6)
            points.Draw("P")
            legend.AddEntry(points, "all runs", "p")
        for edge in (-common.HALF_WINDOW["centroid"], common.HALF_WINDOW["centroid"]):
            line = common.keep(ROOT.TLine(edge, frame.GetYaxis().GetXmin(), edge, frame.GetYaxis().GetXmax()))
            line.Draw()
        legend.Draw()
        if fit is not None:
            fit["function"].SetLineColor(ROOT.kRed)
            fit["function"].SetLineWidth(2)
            common.keep(fit["function"]).Draw("same")
            lines = [f"c = {fit['coefficients'][2]:.1f} #pm {fit['errors'][2]:.1f} ADC/cr^{{2}}   "
                     f"c/a = {fit['rel']:.2f} #pm {fit['err_rel']:.2f} %/cr^{{2}}   "
                     f"#chi^{{2}}/ndf {fit['chi2']:.1f} / {fit['ndf']}"]
            if raw_fit is not None and len(runs) > 1:
                lines.append(f"without normalisation: c/a = {raw_fit['rel']:.2f} %/cr^{{2}}   "
                             f"#chi^{{2}}/ndf {raw_fit['chi2']:.1f} / {raw_fit['ndf']}")
            common.keep(common.text_box(lines, 0.13, 0.89 - 0.045 * len(lines), 0.89, 0.89, 0.026)).Draw()
            rows.append(dict(resistance=resistance, energy=energy, coord=name, nrun=len(runs),
                             nev=int(in_profile.sum()),
                             a=fit["coefficients"][0], b=fit["coefficients"][1], c=fit["coefficients"][2],
                             err_a=fit["errors"][0], err_b=fit["errors"][1], err_c=fit["errors"][2],
                             rel_pct=fit["rel"], err_rel_pct=fit["err_rel"], chi2=fit["chi2"], ndf=fit["ndf"],
                             rel_pct_raw=raw_fit["rel"] if raw_fit else np.nan,
                             chi2_raw=raw_fit["chi2"] if raw_fit else np.nan,
                             ndf_raw=raw_fit["ndf"] if raw_fit else np.nan))
    common.save_canvas(canvas, os.path.join(outdir, "profiles", f"profilo_{energy}GeV_{resistance}ohm.png"))
    return rows


# ------------------------------------------------------------------ maps, drift_dcb_all style
def centroid_maps(events, base_box, fits, pooled, resistance, energy, outdir):
    """Occupancy, <A> map, run centroids and the two pooled profiles with the parabola
    p1 + p2 (x - p0)^2. Returns the parabola rows."""
    threshold = map_threshold(energy, resistance)
    run, amplitude, u_values, v_values = events["run"], events["amplitude"], events["u"], events["v"]
    above = base_box & (amplitude > threshold)
    canvas = ROOT.TCanvas(f"maps_{resistance}_{energy}", "", 1800, 1100)
    canvas.Divide(3, 2)

    pad = canvas.cd(1)
    pad.SetRightMargin(0.16)
    occupancy = common.keep(ROOT.TH2D(common.unique_name("occupancy"),
                                      f"{resistance} #Omega, {energy} GeV: occupancy in |pos_eta - 18|, |pos_phi - 6| < {MAP_HALF}, N = {int(base_box.sum())};pos_eta;pos_phi",
                                      MAP_NBINS, common.ETA_CENTRE - MAP_HALF, common.ETA_CENTRE + MAP_HALF,
                                      MAP_NBINS, common.PHI_CENTRE - MAP_HALF, common.PHI_CENTRE + MAP_HALF))
    occupancy.FillN(int(base_box.sum()), np.ascontiguousarray(events["pos_eta"][base_box]),
                    np.ascontiguousarray(events["pos_phi"][base_box]), np.ones(int(base_box.sum())))
    occupancy.Draw("COLZ")
    centroid = common.keep(ROOT.TMarker(events["pos_eta"][base_box].mean(), events["pos_phi"][base_box].mean(), 34))
    centroid.SetMarkerColor(ROOT.kRed)
    centroid.SetMarkerSize(2)
    centroid.Draw()

    pad = canvas.cd(2)
    pad.SetRightMargin(0.16)
    counts, x_edges, y_edges = np.histogram2d(events["pos_eta"][above], events["pos_phi"][above],
                                              bins=MAP_NBINS,
                                              range=[[common.ETA_CENTRE - MAP_HALF, common.ETA_CENTRE + MAP_HALF],
                                                     [common.PHI_CENTRE - MAP_HALF, common.PHI_CENTRE + MAP_HALF]])
    sums, _, _ = np.histogram2d(events["pos_eta"][above], events["pos_phi"][above], bins=[x_edges, y_edges],
                                weights=amplitude[above])
    inclusive = amplitude[above].mean() if above.any() else 1.
    mean_map = common.keep(ROOT.TH2D(common.unique_name("mean_map"),
                                     f"#LTA#GT, A > {threshold:.0f} ADC, N_{{bin}} #geq {MAP_MIN_PER_BIN};pos_eta;pos_phi",
                                     MAP_NBINS, x_edges[0], x_edges[-1], MAP_NBINS, y_edges[0], y_edges[-1]))
    for x_index in range(MAP_NBINS):
        for y_index in range(MAP_NBINS):
            if counts[x_index, y_index] >= MAP_MIN_PER_BIN:
                mean_map.SetBinContent(x_index + 1, y_index + 1, sums[x_index, y_index] / counts[x_index, y_index])
    mean_map.SetMinimum(inclusive * (1 - PALETTE_HALF))      # empty bins (0) stay white
    mean_map.SetMaximum(inclusive * (1 + PALETTE_HALF))
    mean_map.Draw("COLZ")

    pad = canvas.cd(3)
    runs = sorted(fits)
    x_index = np.arange(len(runs), dtype=float)
    per_run_masks = [above & (run == this_run) for this_run in runs]
    eta_mean = np.array([events["pos_eta"][mask].mean() if mask.any() else np.nan for mask in per_run_masks])
    phi_mean = np.array([events["pos_phi"][mask].mean() if mask.any() else np.nan for mask in per_run_masks])
    eta_err = np.array([events["pos_eta"][mask].std() / np.sqrt(max(mask.sum(), 1)) for mask in per_run_masks])
    phi_err = np.array([events["pos_phi"][mask].std() / np.sqrt(max(mask.sum(), 1)) for mask in per_run_masks])
    global_eta = events["pos_eta"][above].mean() if above.any() else 0.
    global_phi = events["pos_phi"][above].mean() if above.any() else 0.
    if len(runs):
        span = float(np.nanmax(np.abs(np.concatenate([eta_mean - global_eta, phi_mean - global_phi,
                                                      eta_err, phi_err, [1e-3]]))))
        frame = common.run_axis_frame(runs, -1.6 * span, 1.6 * span,
                                      "run centroid - global centroid;run;[crystal units]")
        frame.Draw()
        common.keep(ROOT.TLine(-0.5, 0., len(runs) - 0.5, 0.)).Draw()
        shift_eta = common.keep(common.make_graph(x_index, eta_mean - global_eta, eta_err))
        shift_phi = common.keep(common.make_graph(x_index + 0.1, phi_mean - global_phi, phi_err))
        for graph, marker, colour in ((shift_eta, 20, ROOT.kAzure + 2), (shift_phi, 21, ROOT.kOrange + 7)):
            graph.SetMarkerStyle(marker)
            graph.SetMarkerColor(colour)
            graph.SetLineColor(colour)
            graph.Draw("P")
        legend = common.keep(ROOT.TLegend(0.6, 0.78, 0.88, 0.88))
        legend.AddEntry(shift_eta, "pos_eta", "p")
        legend.AddEntry(shift_phi, "pos_phi", "p")
        legend.Draw()

    parabola_rows = []
    response_window = ((pooled["peak"] - N_SIGMA_WINDOW * pooled["sigma"], pooled["peak"] + N_SIGMA_WINDOW * pooled["sigma"])
                       if pooled else (0.8 * common.SCALE[resistance] * energy, 1.2 * common.SCALE[resistance] * energy))
    for pad_index, (name, own, other, label, other_cut) in enumerate(COORDINATES):
        pad = canvas.cd(4 + pad_index)
        pad.SetLeftMargin(0.13)
        in_profile = (above & (np.abs(events[other]) < common.HALF_WINDOW["centroid"])
                      & (amplitude > response_window[0]) & (amplitude < response_window[1]))
        centres, means, errors, _n = binned_profile(events[own][in_profile], amplitude[in_profile])
        if len(centres) < 4:
            continue
        points = common.keep(common.make_graph(centres, means, errors))
        points.SetTitle(f"{other_cut},  {response_window[0]:.0f} < A < {response_window[1]:.0f} ADC;"
                        f"{label} [crystal units];#LTA#GT [ADC]")
        points.SetMarkerStyle(20)
        points.SetMarkerSize(0.6)
        points.Draw("AP")
        points.GetXaxis().SetLimits(-PROFILE_HALF, PROFILE_HALF)
        padding = 0.12 * (means.max() - means.min()) + 2 * np.median(errors)
        points.GetYaxis().SetRangeUser(means.min() - padding, means.max() + padding)
        for edge in (-common.HALF_WINDOW["centroid"], common.HALF_WINDOW["centroid"]):
            common.keep(ROOT.TLine(edge, means.min() - padding, edge, means.max() + padding)).Draw()
        inside = np.abs(centres) <= FIT_HALF
        if inside.sum() >= 5:
            fit_graph = common.keep(common.make_graph(centres[inside], means[inside], errors[inside]))
            parabola = common.keep(ROOT.TF1(common.unique_name("vertex_parabola"), "[1] + [2]*(x-[0])^2",
                                            centres[inside].min(), centres[inside].max()))
            parabola.SetParNames("p0", "p1", "p2")
            parabola.SetParameters(0., float(means[inside].max()), -100.)
            result = common.fit_graph(fit_graph, parabola, lower_limit_zero=False)
            parabola.SetLineColor(ROOT.kRed)
            parabola.Draw("same")
            values, errs = result["values"], result["errors"]
            ndf = int(inside.sum()) - 3
            common.keep(common.text_box(
                [f"fit over |x| < {FIT_HALF}   #chi^{{2}}/ndf {result['chi2']:.1f} / {ndf}",
                 f"p0 (vertex) {values[0]:+.4f} #pm {errs[0]:.4f}   p1 (max) {values[1]:.1f} #pm {errs[1]:.1f}",
                 f"p2 (curv.) {values[2]:.1f} #pm {errs[2]:.1f}"],
                0.14, 0.12, 0.9, 0.27, 0.028)).Draw()
            parabola_rows.append(dict(resistance=resistance, energy=energy, coord=name, n_bin=int(inside.sum()),
                                      p0_vertice=values[0], err_p0=errs[0], p1_massimo=values[1], err_p1=errs[1],
                                      p2_curvatura=values[2], err_p2=errs[2], chi2=result["chi2"], ndf=ndf))
    common.save_canvas(canvas, os.path.join(outdir, "maps", f"centroide2D_{energy}GeV_{resistance}ohm.png"))
    return parabola_rows


# ------------------------------------------------------------------ summaries
def curvature_summary(profile_rows, resistances, outdir):
    canvas = ROOT.TCanvas("curvature", "", 1300, 500)
    canvas.Divide(2, 1)
    for pad_index, (name, _own, _other, label, _cut) in enumerate(COORDINATES):
        pad = canvas.cd(pad_index + 1)
        pad.SetLogx()
        pad.SetGrid()
        legend = common.keep(ROOT.TLegend(0.15, 0.7, 0.4, 0.88))
        first = True
        for resistance in resistances:
            rows = sorted([row for row in profile_rows if row["resistance"] == resistance and row["coord"] == name],
                          key=lambda row: row["energy"])
            if not rows:
                continue
            graph = common.keep(common.make_graph([row["energy"] for row in rows], [-row["rel_pct"] for row in rows],
                                                  [row["err_rel_pct"] for row in rows]))
            graph.SetMarkerStyle(common.MARKER[resistance])
            graph.SetMarkerColor(common.COLOUR[resistance])
            graph.SetLineColor(common.COLOUR[resistance])
            graph.SetTitle(f"{label};E_{{nom}} [GeV];-c/a  [% / crystal^{{2}}]")
            graph.Draw("APL" if first else "PL")
            if first:
                graph.GetXaxis().SetLimits(15, 320)
                graph.GetYaxis().SetRangeUser(0, 1.3 * max(-row["rel_pct"] for row in profile_rows if row["coord"] == name))
            legend.AddEntry(graph, f"{resistance} #Omega", "pl")
            first = False
        legend.Draw()
    common.save_canvas(canvas, os.path.join(outdir, "profiles", "curvatura_vs_energia.png"))

    for resistance in resistances:
        rows = [row for row in profile_rows if row["resistance"] == resistance]
        if not rows:
            continue
        energies = sorted({row["energy"] for row in rows})
        canvas = ROOT.TCanvas(f"parabolas_{resistance}", "", 1500, 560)
        canvas.Divide(2, 1)
        for pad_index, (name, _own, _other, label, other_cut) in enumerate(COORDINATES):
            pad = canvas.cd(pad_index + 1)
            pad.SetGrid()
            frame = common.keep(ROOT.TH2F(common.unique_name("frame"), f"{other_cut};{label} [crystal units];(a + bx + cx^{{2}}) / a",
                                          10, -FIT_HALF, FIT_HALF, 10, 0.96, 1.005))
            frame.Draw()
            legend = common.keep(ROOT.TLegend(0.14, 0.12, 0.34, 0.12 + 0.03 * len(energies)))
            legend.SetTextSize(0.025)
            for energy_index, energy in enumerate(energies):
                match = [row for row in rows if row["coord"] == name and row["energy"] == energy]
                if not match:
                    continue
                row = match[0]
                curve = common.keep(ROOT.TF1(common.unique_name("norm_parabola"),
                                             f"({row['a']} + {row['b']}*x + {row['c']}*x*x) / {row['a']}", -FIT_HALF, FIT_HALF))
                curve.SetLineColor(ROOT.TColor.GetPalette()[int(240 * energy_index / max(len(energies) - 1, 1))])
                curve.Draw("same")
                legend.AddEntry(curve, f"{energy} GeV", "l")
            legend.Draw()
        common.save_canvas(canvas, os.path.join(outdir, "profiles", f"profili_all_energies_{resistance}ohm.png"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", required=True)
    parser.add_argument("--workdir", required=True, help="the --outdir of stage 1 (centroid pass)")
    parser.add_argument("--resistances", nargs="+", type=int, default=[340, 400, 500])
    parser.add_argument("--exclude-runs", nargs="*", type=int, default=[])
    runsets.add_argument(parser)
    args = parser.parse_args()
    common.style()

    fit_rows = common.read_csv(common.require(os.path.join(args.workdir, "01_dcb_per_run.csv")))
    windows = {(row["resistance"], row["energy"]): row
               for row in common.read_csv(common.require(os.path.join(args.workdir, "01_windows.csv")))}
    if fit_rows and fit_rows[0]["selection"] != "centroid":
        parser.error("stage 3 belongs to the centroid pass")
    amplitude = fit_rows[0]["amplitude"] if fit_rows else "a3x3"
    dropped, kept_only = runsets.resolve(args.runset, args.exclude_runs)
    excluded = {key for key in windows if windows[key]["skipped"]}

    moment_rows, profile_rows, parabola_rows = [], [], []
    for resistance, energy, path in common.resistance_energy_pairs(args.base, args.resistances, excluded):
        if (resistance, energy) not in windows:
            continue
        fits = per_run_fits(fit_rows, resistance, energy)
        if not fits:
            print(f"[{resistance} ohm {energy:>4} GeV] no per-run fit, skipped")
            continue
        print(f"[{resistance} ohm {energy:>4} GeV] {len(fits)} runs", flush=True)
        events = common.read_events(path, amplitude)
        half = windows[(resistance, energy)]["half"]
        base = (events["A_tot"] > common.A_TOT_MIN) & common.runset_mask(events["run"], dropped, kept_only)
        selected = base & (np.abs(events["u"]) <= half) & (np.abs(events["v"]) <= half)
        moment_rows += moments_rows(events, selected, fits, resistance, energy)
        profile_rows += normalised_profiles(events, base, fits, resistance, energy, args.workdir)
        pooled = [row for row in fit_rows if row["resistance"] == resistance and row["energy"] == energy
                  and row["variation"] == "nominal" and row["run"] == 0]
        in_box = (common.runset_mask(events["run"], dropped, kept_only)
                  & (np.abs(events["u"]) < MAP_HALF) & (np.abs(events["v"]) < MAP_HALF))
        if in_box.sum() > 500:
            parabola_rows += centroid_maps(events, in_box, fits, pooled[0] if pooled else None,
                                           resistance, energy, args.workdir)

    common.write_csv(os.path.join(args.workdir, "03_moments.csv"), moment_rows, MOMENT_COLUMNS)
    common.write_csv(os.path.join(args.workdir, "03_profiles.csv"), profile_rows, PROFILE_COLUMNS)
    common.write_csv(os.path.join(args.workdir, "03_parabola_centroide.csv"), parabola_rows, PARABOLA_COLUMNS)
    if profile_rows:
        curvature_summary(profile_rows, args.resistances, args.workdir)


if __name__ == "__main__":
    main()
