#!/usr/bin/env python3
"""
Stage 9 -- the resolution points and the size of every term, and save.

Reads 08_systematics.csv and draws, for the selection of the pass, one column per
resistance: on top sigma/mu as measured and sigma/E after the subtractions with its
total error bar; below, on a log scale, every term that entered (BES, synchrotron,
drift, and the recipe's own systematics), with markers for the points where the drift
is undefined (one run) or exactly zero (chi2/ndf <= 1). Pooled points and points with a
hand-set hodoscope window are marked as such. No fit here: stage 10 does the fits.

Writes 09_resolution_points.csv (energy_true, sigma_over_E, its error, and in_fit, the
flag stage 10 uses: 0 for the energies listed in --nofit-energies) and
09_resolution_terms.png. With --bes nominal (codiceA only) the central value has the
nominal BES (BES_formula) subtracted instead of the conservative one, the error bar has
no BES systematic (err_total_nominal_bes), and the files carry the suffix _nominal_bes.
"""

import argparse
import os

import numpy as np
import ROOT

import common

COLUMNS = ("resistance", "energy", "energy_true", "selection", "recipe", "window", "n_run", "pooled",
           "sigma_raw", "sigma_over_E", "err", "in_fit")
TERM_STYLE = {           # column: (legend, marker, colour)
    "bes": ("BES", 27, ROOT.kOrange + 7), "sync": ("synchrotron", 26, ROOT.kViolet),
    "drift": ("drift", 21, ROOT.kAzure + 2), "vtx_syst": ("vertex syst", 25, ROOT.kTeal + 3),
    "bes_syst": ("BES syst", 28, ROOT.kMagenta + 2), "syn_syst": ("sync syst", 30, ROOT.kGray + 2),
    "unif_syst": ("map syst", 34, ROOT.kGreen + 2), "stat": ("stat", 24, ROOT.kBlack)}
TERMS_BY_RECIPE = {"codiceA": ("bes", "sync", "drift", "stat", "vtx_syst", "bes_syst", "syn_syst"),
                   "uniforme": ("bes", "sync", "drift", "stat", "unif_syst")}


def graph_of(rows, key, error_key=None, positive_only=False):
    x_values = np.array([row["energy_true"] for row in rows], float)
    y_values = np.array([row[key] for row in rows], float)
    errors = np.array([row[error_key] for row in rows], float) if error_key else None
    keep = np.isfinite(y_values) & ((y_values > 0) if positive_only else True)
    if not keep.any():
        return None
    return common.make_graph(x_values[keep], y_values[keep], errors[keep] if errors is not None else None)


def styled(graph, marker, colour, size=1.1):
    graph.SetMarkerStyle(marker)
    graph.SetMarkerColor(colour)
    graph.SetLineColor(colour)
    graph.SetMarkerSize(size)
    return graph


def draw_top(pad, rows, resistance, selection, central="sigma_corr", bes_label="BES", error_column="err_total"):
    pad.SetGrid()
    pad.SetLeftMargin(0.13)
    raw = common.keep(styled(graph_of(rows, "sigma_raw"), 20, ROOT.kGray + 2))
    raw.SetTitle(f"{resistance} #Omega, cut on the {selection};E_{{true}} [GeV];#sigma/E  [%]")
    raw.Draw("APL")
    raw.GetXaxis().SetLimits(0.85 * min(row["energy_true"] for row in rows), 1.15 * max(row["energy_true"] for row in rows))
    raw.GetYaxis().SetRangeUser(0, 1.08 * max(row["sigma_raw"] for row in rows))
    legend = common.keep(ROOT.TLegend(0.45, 0.6, 0.89, 0.88))
    legend.SetTextSize(0.03)
    legend.AddEntry(raw, "#sigma/#mu", "pl")
    corrected = graph_of(rows, central, error_column, positive_only=True)
    if corrected is not None:
        common.keep(styled(corrected, 22, ROOT.kRed + 1, 1.4)).Draw("PL")
        legend.AddEntry(corrected, f"- {bes_label} - synchrotron" if rows[0]["recipe"] == "uniforme"
                        else f"- {bes_label} - synchrotron, err with systs", "pl")
    pooled = [row for row in rows if row["pooled"] and np.isfinite(row[central])]
    if pooled:
        marker = common.keep(styled(graph_of(pooled, central), 25, ROOT.kViolet, 2.4))
        marker.Draw("P")
        legend.AddEntry(marker, "one pooled fit (no run has enough events)", "p")
    hand_set = [row for row in rows if str(row["window"]).endswith("-ecal_prof") and np.isfinite(row[central])]
    if hand_set:
        marker = common.keep(styled(graph_of(hand_set, central), 24, ROOT.kGray + 3, 2.4))
        marker.Draw("P")
        legend.AddEntry(marker, "no parabola: hand-set hodoscope window", "p")
    legend.Draw()


def draw_terms(pad, rows, recipe, bes_column="bes"):
    pad.SetGrid()
    pad.SetLogy()
    pad.SetLeftMargin(0.13)
    energies = [row["energy_true"] for row in rows]
    values = [row["sigma_raw"] for row in rows]
    for term in TERMS_BY_RECIPE[recipe]:
        values += [row[term] for row in rows if np.isfinite(row[term]) and row[term] > 0]
    low, high = 0.12 * min(values), 2.5 * max(values)
    frame = common.keep(ROOT.TH2F(common.unique_name("terms_frame"), ";E_{true} [GeV];size of each term  [%]",
                                  10, 0.85 * min(energies), 1.15 * max(energies), 10, low, high))
    frame.Draw()
    legend = common.keep(ROOT.TLegend(0.14, 0.12, 0.55, 0.12 + 0.032 * (len(TERMS_BY_RECIPE[recipe]) + 3)))
    legend.SetTextSize(0.028)
    raw = common.keep(styled(graph_of(rows, "sigma_raw"), 20, ROOT.kGray + 2))
    raw.Draw("PL")
    legend.AddEntry(raw, "#sigma/#mu", "pl")
    for term in TERMS_BY_RECIPE[recipe]:
        label, marker, colour = TERM_STYLE[term]
        graph = graph_of(rows, bes_column if term == "bes" else term, positive_only=True)
        if graph is None:
            continue
        common.keep(styled(graph, marker, colour)).Draw("PL")
        legend.AddEntry(graph, label, "pl")
    single = [row for row in rows if row["drift"] <= 0 and row["n_run"] <= 1]
    compatible = [row for row in rows if row["drift"] <= 0 and row["n_run"] > 1]
    for subset, marker, label in ((single, 5, "drift n/a (1 run)"), (compatible, 25, "drift = 0 (#chi^{2}/ndf #leq 1)")):
        if subset:
            graph = common.keep(styled(common.make_graph([row["energy_true"] for row in subset],
                                                         [low * 2.] * len(subset)), marker, ROOT.kAzure + 2))
            graph.Draw("P")
            legend.AddEntry(graph, label, "p")
    legend.Draw()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--nofit-energies", nargs="*", type=int, default=[],
                        help="energies drawn but left out of the N/S/C fit of stage 10")
    parser.add_argument("--bes", choices=("cons", "nominal"), default="cons",
                        help="codiceA only: which BES is subtracted in the central value")
    args = parser.parse_args()
    common.style()

    rows = common.read_csv(common.require(os.path.join(args.workdir, "08_systematics.csv")))
    if not rows:
        raise SystemExit("08_systematics.csv is empty")
    selection, recipe = rows[0]["selection"], rows[0]["recipe"]
    if args.bes == "nominal" and recipe != "codiceA":
        parser.error("--bes nominal exists only in the codiceA recipe")
    central = "sigma_corr_nominal_bes" if args.bes == "nominal" else "sigma_corr"
    error_column = "err_total_nominal_bes" if args.bes == "nominal" else "err_total"
    bes_column = "bes_nom" if args.bes == "nominal" else "bes"
    suffix = "_nominal_bes" if args.bes == "nominal" else ""
    points = [dict(resistance=row["resistance"], energy=row["energy"], energy_true=row["energy_true"],
                   selection=selection, recipe=recipe, window=row["window"], n_run=row["n_run"],
                   pooled=row["pooled"], sigma_raw=row["sigma_raw"], sigma_over_E=row[central],
                   err=row[error_column], in_fit=int(row["energy"] not in args.nofit_energies))
              for row in rows]
    common.write_csv(os.path.join(args.workdir, f"09_resolution_points{suffix}.csv"), points, COLUMNS)

    resistances = sorted({row["resistance"] for row in rows})
    canvas = ROOT.TCanvas("resolution_terms", "", 640 * len(resistances), 1000)
    canvas.Divide(len(resistances), 2)
    for column, resistance in enumerate(resistances):
        subset = sorted([row for row in rows if row["resistance"] == resistance], key=lambda row: row["energy"])
        draw_top(canvas.cd(column + 1), subset, resistance, selection, central,
                 "BES (nominal)" if args.bes == "nominal" else "BES", error_column)
        draw_terms(canvas.cd(len(resistances) + column + 1), subset, recipe, bes_column)
    common.save_canvas(canvas, os.path.join(args.workdir, f"09_resolution_terms{suffix}.png"))


if __name__ == "__main__":
    main()
