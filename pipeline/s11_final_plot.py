#!/usr/bin/env python3
"""
Stage 11 -- the final figure, everything on one canvas (ROOT).

One column per resistance. On top: sigma/mu as measured, sigma/E after the
subtractions with its total error, the points with a hand-set hodoscope window or a
pooled fit marked, and the N/S/C fit curves of stage 10 with their parameter boxes
(for the codiceA recipe the fit with S held at the 340 ohm value and the one with S and
C both held; for the uniforme recipe the per-resistance fit and the common one). Below,
on a log scale, the size of every term that entered. Nothing is fitted here: the curves
are drawn from 10_resolution_fits.csv.

The legends sit where the data are not: the top-right corner above the falling curve,
and, for the terms, a strip under the axis.

Reads 08_systematics.csv, 09_resolution_points[_nominal_bes].csv and
10_resolution_fits[_nominal_bes].csv; writes 11_resolution[_nominal_bes].png (and the
canvas in 11_resolution[_nominal_bes].root). s11_final_plot_matplotlib.py draws the same
figure with matplotlib.
"""

import argparse
import os

import numpy as np
import ROOT

import common
from s09_resolution_plots import TERMS_BY_RECIPE, TERM_STYLE, graph_of, styled

RESISTANCES = (340, 400, 500)
TERMS_FLOOR = 3e-4          # the terms axis does not go below this (the sync syst does)
# one line style per fit variant: (legend label, line style, colour)
FIT_STYLE = {("indep", "S"): ("fit, S from 340 #Omega", 2, ROOT.kViolet),
             ("indep", "SC"): ("fit, S and C from 340 #Omega", 3, ROOT.kGreen + 2),
             ("indep", ""): ("fit N/E #oplus S/#sqrt{E} #oplus C", 2, ROOT.kViolet),
             ("common", ""): ("fit, S and C common", 3, ROOT.kGreen + 2)}


def parameter_lines(row):
    def one(name, value, error, fixed, unit, digits):
        return f"{name}  {value:{digits}}" + (f" {unit} (fixed)" if fixed else f" #pm {error:{digits}} {unit}")
    return [one("N", row["N_MeV"], row["err_N_MeV"], row["N_fixed"], "MeV", "5.0f"),
            one("S", row["S_pct"], row["err_S"], row["S_fixed"], "%", "6.3f"),
            one("C", row["C_pct"], row["err_C"], row["C_fixed"], "%", "6.3f"),
            f"#chi^{{2}}/ndf  {row['chi2']:.1f} / {row['ndf']}"]


def draw_top(pad, points, terms, fits, resistance, selection, central_label):
    pad.SetGrid()
    pad.SetLeftMargin(0.13)
    pad.SetTopMargin(0.08)
    energies = [row["energy_true"] for row in points]
    top = 1.9 * max(row["sigma_raw"] for row in points)          # room for legend and boxes
    frame = common.keep(ROOT.TH2F(common.unique_name("top_frame"),
                                  f"{resistance} #Omega, cut on the {selection}, A_{{tot}} > {common.A_TOT_MIN:.0f} ADC"
                                  f";E_{{true}} [GeV];#sigma/E  [%]",
                                  10, 0.85 * min(energies), 1.15 * max(energies), 10, 0., top))
    frame.Draw()
    legend = common.keep(ROOT.TLegend(0.45, 0.62, 0.89, 0.9))
    legend.SetTextSize(0.03)
    legend.SetBorderSize(1)
    raw = common.keep(styled(graph_of(points, "sigma_raw"), 20, ROOT.kGray + 2))
    raw.Draw("PL")
    legend.AddEntry(raw, "#sigma/#mu", "pl")
    corrected = common.keep(styled(graph_of(points, "sigma_over_E", "err", positive_only=True), 22, ROOT.kRed + 1, 1.4))
    corrected.Draw("PL")
    legend.AddEntry(corrected, central_label, "pl")
    hand_set = [row for row in points if str(row["window"]).endswith("-ecal_prof")]
    if hand_set:
        marker = common.keep(styled(graph_of(hand_set, "sigma_over_E"), 24, ROOT.kGray + 3, 2.4))
        marker.Draw("P")
        legend.AddEntry(marker, "no parabola: hand-set hodoscope window", "p")
    pooled = [row for row in points if row["pooled"]]
    if pooled:
        marker = common.keep(styled(graph_of(pooled, "sigma_over_E"), 25, ROOT.kViolet, 2.4))
        marker.Draw("P")
        legend.AddEntry(marker, "one pooled fit (no run has enough events)", "p")
    left_out = [row for row in points if not row["in_fit"]]
    if left_out:
        marker = common.keep(styled(graph_of(left_out, "sigma_over_E"), 26, ROOT.kGray + 2, 2.4))
        marker.Draw("P")
        legend.AddEntry(marker, "not in the fit", "p")
    box_top = 0.6
    for fit in fits:
        label, style, colour = FIT_STYLE[(fit["mode"], str(fit["fixed_from_340"]))]
        curve = common.resolution_function(common.unique_name("curve"), 0.9 * min(energies), 1.05 * max(energies))
        curve.SetParameters(fit["N_MeV"] / 1000., fit["S_pct"], fit["C_pct"])
        curve.SetLineColor(colour)
        curve.SetLineStyle(style)
        curve.SetLineWidth(3)
        common.keep(curve).Draw("same")
        legend.AddEntry(curve, label, "l")
        common.keep(common.text_box(parameter_lines(fit), 0.5, box_top - 0.17, 0.89, box_top, 0.028, colour)).Draw()
        box_top -= 0.19
    legend.Draw()


def draw_terms(pad, terms, recipe, bes_column):
    pad.SetGrid()
    pad.SetLogy()
    pad.SetLeftMargin(0.13)
    pad.SetBottomMargin(0.36)                      # the legend lives under the axis
    energies = [row["energy_true"] for row in terms]
    values = [row["sigma_raw"] for row in terms]
    for term in TERMS_BY_RECIPE[recipe]:
        column = bes_column if term == "bes" else term
        values += [row[column] for row in terms if np.isfinite(row[column]) and row[column] > 0]
    low, high = max(0.25 * min(values), TERMS_FLOOR), 2.5 * max(values)
    frame = common.keep(ROOT.TH2F(common.unique_name("terms_frame"), ";E_{true} [GeV];size of each term  [%]",
                                  10, 0.85 * min(energies), 1.15 * max(energies), 10, low, high))
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.Draw()
    legend = common.keep(ROOT.TLegend(0.13, 0.01, 0.9, 0.22))
    legend.SetNColumns(3)
    legend.SetTextSize(0.028)
    legend.SetBorderSize(0)
    raw = common.keep(styled(graph_of(terms, "sigma_raw"), 20, ROOT.kGray + 2))
    raw.Draw("PL")
    legend.AddEntry(raw, "#sigma/#mu", "pl")
    for term in TERMS_BY_RECIPE[recipe]:
        label, marker, colour = TERM_STYLE[term]
        graph = graph_of(terms, bes_column if term == "bes" else term, positive_only=True)
        if graph is None:
            continue
        common.keep(styled(graph, marker, colour)).Draw("PL")
        legend.AddEntry(graph, label, "pl")
    single = [row for row in terms if row["drift"] <= 0 and row["n_run"] <= 1]
    compatible = [row for row in terms if row["drift"] <= 0 and row["n_run"] > 1]
    for subset, marker, label in ((single, 5, "drift n/a (1 run)"), (compatible, 25, "drift = 0 (#chi^{2}/ndf #leq 1)")):
        if subset:
            graph = common.keep(styled(common.make_graph([row["energy_true"] for row in subset],
                                                         [low * 1.6] * len(subset)), marker, ROOT.kAzure + 2))
            graph.Draw("P")
            legend.AddEntry(graph, label, "p")
    legend.Draw()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--bes", choices=("cons", "nominal"), default="cons")
    args = parser.parse_args()
    common.style()
    suffix = "_nominal_bes" if args.bes == "nominal" else ""

    terms_all = common.read_csv(common.require(os.path.join(args.workdir, "08_systematics.csv")))
    points_all = common.read_csv(common.require(os.path.join(args.workdir, f"09_resolution_points{suffix}.csv")))
    fits_all = [row for row in common.read_csv(common.require(os.path.join(args.workdir, f"10_resolution_fits{suffix}.csv")))
                if row["variant"] == "nominal"]
    selection, recipe = points_all[0]["selection"], points_all[0]["recipe"]
    resistances = [resistance for resistance in RESISTANCES if any(row["resistance"] == resistance for row in points_all)]
    bes_column = "bes_nom" if args.bes == "nominal" else "bes"
    central_label = ("- BES - synchrotron" if recipe == "uniforme"
                     else f"- BES ({'nominal' if args.bes == 'nominal' else 'conservative'}) - synchrotron")

    canvas = ROOT.TCanvas("final", "", 700 * len(resistances), 1150)
    canvas.Divide(len(resistances), 2)
    for column, resistance in enumerate(resistances):
        points = sorted([row for row in points_all if row["resistance"] == resistance], key=lambda row: row["energy"])
        terms = sorted([row for row in terms_all if row["resistance"] == resistance], key=lambda row: row["energy"])
        fits = [row for row in fits_all if row["resistance"] == resistance]
        draw_top(canvas.cd(column + 1), points, terms, fits, resistance, selection, central_label)
        draw_terms(canvas.cd(len(resistances) + column + 1), terms, recipe, bes_column)
    common.save_canvas(canvas, os.path.join(args.workdir, f"11_resolution{suffix}.png"))
    output = ROOT.TFile(os.path.join(args.workdir, f"11_resolution{suffix}.root"), "RECREATE")
    canvas.Write("final")
    output.Close()


if __name__ == "__main__":
    main()
