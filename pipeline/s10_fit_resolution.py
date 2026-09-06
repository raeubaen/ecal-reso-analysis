#!/usr/bin/env python3
"""
Stage 10 -- the fits of sigma/E = N/E (+) S/sqrt(E) (+) C, and save.

Reads 09_resolution_points.csv and fits with ROOT (TGraphErrors + TF1, MIGRAD + HESSE)
following the recipe of the selection:

  codiceA (hodoscope, resolution_hodo.py)
      N, S, C >= 0 per resistance, seeds N = 0.3 GeV, S = 3 %, C = 0.3 %. 340 ohm is
      fitted first; for the other resistances two fits are done, one with S held at the
      340 ohm value (fixed_from_340 = S, the flat script, C free) and one with S and C
      both held at the 340 ohm values (fixed_from_340 = SC, only N free). (The boxes of the slides of 10 September 2026 print "C (FIXED)" for 400
      and 500 ohm, but that label is applied to both lines whenever R != 340 in the flat
      script: C was free, as its different values show.) Two variants when
      --nofit-energies was given to stage 9: with every point, and with those energies
      drawn but left out. --bes nominal reads the points of stage 9 --bes nominal.
  uniforme (centroid, resolution_final_uniforme.py)
      N, S, C >= 0 per resistance, seeds N = 0.3, S = 5 (the 340 ohm S at 500 ohm),
      C = 0.3; at 500 ohm C is held at the seed. Then the simultaneous fit with S and C
      common to the three resistances and one N each, done as in root/fit_resolution.C:
      one combined TGraphErrors with abscissa E + 10000 * resistance index and a
      5-parameter TF1 that decodes it.

Writes 10_resolution_fits.csv (N in MeV, S and C in %, chi2, ndf, what was held fixed),
10_resolution_fits.root (graphs, functions, canvases) and the canvases
10_resolution_fit[_allpoints].png and 10_resolution_common.png.
"""

import argparse
import os

import numpy as np
import ROOT

import common

COLUMNS = ("selection", "recipe", "variant", "mode", "fixed_from_340", "resistance", "n_points",
           "N_MeV", "err_N_MeV", "N_fixed", "S_pct", "err_S", "S_fixed", "C_pct", "err_C", "C_fixed",
           "chi2", "ndf")
FIXED_SETS = {"codiceA": ("S", "SC"), "uniforme": ("",)}
RESISTANCES = (340, 400, 500)
X_OFFSET = 10000.       # resistance index encoded in the abscissa of the combined graph
C_FIXED_UNIFORME = 0.3  # C at 500 ohm in resolution_final_uniforme.py (held at its seed)


def points_of(rows, resistance, fit_only=True):
    subset = sorted([row for row in rows if row["resistance"] == resistance
                     and np.isfinite(row["sigma_over_E"]) and row["sigma_over_E"] > 0 and row["err"] > 0
                     and (row["in_fit"] or not fit_only)], key=lambda row: row["energy_true"])
    return subset


def fit_one(rows, resistance, recipe, from_340, fixed_from_340):
    """Per-resistance fit; from_340 = {parameter: value} of the 340 ohm fit, fixed_from_340
    the letters of the parameters held at those values for 400 and 500 ohm."""
    if len(rows) < 4:
        return None
    graph = common.make_graph([row["energy_true"] for row in rows], [row["sigma_over_E"] for row in rows],
                              [row["err"] for row in rows])
    low, high = 0.9 * rows[0]["energy_true"], 1.05 * rows[-1]["energy_true"]
    function = common.resolution_function(common.unique_name("reso"), low, high)
    fixed = []
    if recipe == "codiceA":
        function.SetParameters(0.3, 3., 0.3)
        if resistance != 340 and from_340:
            for name in fixed_from_340:
                function.SetParameter({"N": 0, "S": 1, "C": 2}[name], from_340[name])
                fixed.append(name)
        ndf = len(rows) - (2 if resistance == 500 else 3)        # as written in resolution_hodo.py
    else:
        function.SetParameters(0.3, from_340["S"] if (resistance == 500 and from_340) else 5., C_FIXED_UNIFORME)
        if resistance == 500:
            fixed.append("C")
        ndf = len(rows) - (2 if resistance == 500 else 3)
    result = common.fit_graph(graph, function, fixed=fixed)
    return dict(graph=graph, function=function, result=result, ndf=ndf, fixed=fixed)


def fit_common(rows_by_resistance):
    """S and C common, N per resistance, through the combined graph of fit_resolution.C."""
    x_values, y_values, errors = [], [], []
    for index, resistance in enumerate(RESISTANCES):
        for row in rows_by_resistance.get(resistance, []):
            x_values.append(row["energy_true"] + X_OFFSET * index)
            y_values.append(row["sigma_over_E"])
            errors.append(row["err"])
    if len(x_values) < 6:
        return None
    graph = common.make_graph(x_values, y_values, errors)
    formula = ("sqrt(pow(100*((x<10000)*[2] + (x>=10000 && x<20000)*[3] + (x>=20000)*[4])"
               "/(x - 10000*floor(x/10000)), 2) + [0]*[0]/(x - 10000*floor(x/10000)) + [1]*[1])")
    function = ROOT.TF1(common.unique_name("reso_common"), formula, 0., 3 * X_OFFSET)
    function.SetParNames("S", "C", "N340", "N400", "N500")
    function.SetParameters(2.5, 0.35, 0.30, 0.30, 0.21)
    result = common.fit_graph(graph, function)
    return dict(graph=graph, function=function, result=result, ndf=len(x_values) - 5)


def parameter_lines(values, errors, fixed, chi2, ndf):
    return [f"N  {1000 * values[0]:6.0f}" + (" MeV (fixed)" if "N" in fixed else f" #pm {1000 * errors[0]:.0f} MeV"),
            f"S  {values[1]:6.3f}" + (" % (fixed)" if "S" in fixed else f" #pm {errors[1]:.3f} %"),
            f"C  {values[2]:6.3f}" + (" % (fixed)" if "C" in fixed else f" #pm {errors[2]:.3f} %"),
            f"#chi^{{2}}/ndf  {chi2:.1f} / {ndf}"]


def draw_points(pad, rows, all_rows, resistance, selection):
    pad.SetGrid()
    pad.SetLeftMargin(0.13)
    graph = common.keep(common.make_graph([row["energy_true"] for row in all_rows],
                                          [row["sigma_over_E"] for row in all_rows],
                                          [row["err"] for row in all_rows]))
    graph.SetMarkerStyle(22)
    graph.SetMarkerSize(1.4)
    graph.SetMarkerColor(ROOT.kRed + 1)
    graph.SetLineColor(ROOT.kRed + 1)
    graph.SetTitle(f"{resistance} #Omega, cut on the {selection};E_{{true}} [GeV];#sigma/E  [%]")
    graph.Draw("AP")
    graph.GetXaxis().SetLimits(0.85 * all_rows[0]["energy_true"], 1.15 * all_rows[-1]["energy_true"])
    graph.GetYaxis().SetRangeUser(0, 1.3 * max(row["sigma_over_E"] + row["err"] for row in all_rows))
    left_out = [row for row in all_rows if not row["in_fit"]]
    if left_out:
        marker = common.keep(common.make_graph([row["energy_true"] for row in left_out],
                                               [row["sigma_over_E"] for row in left_out]))
        marker.SetMarkerStyle(26)
        marker.SetMarkerSize(2.4)
        marker.SetMarkerColor(ROOT.kGray + 2)
        marker.Draw("P")
    return graph


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--bes", choices=("cons", "nominal"), default="cons",
                        help="codiceA only: fit the points with the nominal BES subtracted")
    args = parser.parse_args()
    common.style()
    suffix = "_nominal_bes" if args.bes == "nominal" else ""

    rows = common.read_csv(common.require(os.path.join(args.workdir, f"09_resolution_points{suffix}.csv")))
    if not rows:
        raise SystemExit("09_resolution_points.csv is empty")
    selection, recipe = rows[0]["selection"], rows[0]["recipe"]
    resistances = [resistance for resistance in RESISTANCES if any(row["resistance"] == resistance for row in rows)]
    output = ROOT.TFile(os.path.join(args.workdir, f"10_resolution_fits{suffix}.root"), "RECREATE")
    table = []

    variants = [("nominal", True)]
    if any(not row["in_fit"] for row in rows):
        variants.append(("allpoints", False))
    for variant, fit_only in variants:
      for fixed_from_340 in FIXED_SETS[recipe]:
        tag = f"{variant}" + (f"_{fixed_from_340}" if fixed_from_340 else "")
        canvas = ROOT.TCanvas(f"fit_{tag}", "", 640 * len(resistances), 560)
        canvas.Divide(len(resistances), 1)
        from_340 = {}
        fitted = {}
        for column, resistance in enumerate(resistances):
            all_rows = points_of(rows, resistance, fit_only=False)
            fit_rows = points_of(rows, resistance, fit_only=fit_only)
            if not all_rows:
                continue
            pad = canvas.cd(column + 1)
            draw_points(pad, fit_rows, all_rows, resistance, selection)
            fit = fit_one(fit_rows, resistance, recipe, from_340, fixed_from_340)
            if fit is None:
                print(f"  {resistance} ohm: {len(fit_rows)} points, no fit")
                continue
            fitted[resistance] = fit
            values, errors = fit["result"]["values"], fit["result"]["errors"]
            if resistance == 340:
                from_340 = dict(N=values[0], S=values[1], C=values[2])
            fit["function"].SetLineColor(ROOT.kViolet)
            fit["function"].SetLineStyle(2)
            fit["function"].SetLineWidth(3)
            fit["function"].Draw("same")
            common.keep(common.text_box(parameter_lines(values, errors, fit["fixed"], fit["result"]["chi2"], fit["ndf"]),
                                        0.45, 0.66, 0.89, 0.88, 0.032, ROOT.kViolet)).Draw()
            table.append(dict(selection=selection, recipe=recipe, variant=variant, mode="indep",
                              fixed_from_340=fixed_from_340, resistance=resistance,
                              n_points=len(fit_rows), N_MeV=1000 * values[0], err_N_MeV=1000 * errors[0],
                              N_fixed=int("N" in fit["fixed"]),
                              S_pct=values[1], err_S=errors[1], S_fixed=int("S" in fit["fixed"]),
                              C_pct=values[2], err_C=errors[2], C_fixed=int("C" in fit["fixed"]),
                              chi2=fit["result"]["chi2"], ndf=fit["ndf"]))
            print(f"  [{tag}] {resistance} ohm: N {1000 * values[0]:.0f} +- {1000 * errors[0]:.0f} MeV"
                  f"{' (fixed)' if 'N' in fit['fixed'] else ''}, "
                  f"S {values[1]:.3f} +- {errors[1]:.3f} %{' (fixed)' if 'S' in fit['fixed'] else ''}, "
                  f"C {values[2]:.3f} +- {errors[2]:.3f} %{' (fixed)' if 'C' in fit['fixed'] else ''}, "
                  f"chi2/ndf {fit['result']['chi2']:.1f}/{fit['ndf']}")
            output.cd()
            fit["graph"].Write(f"gr_{tag}_{resistance}")
            fit["function"].Write(f"f_{tag}_{resistance}_indep")
        output.cd()
        canvas.Write(f"c_{tag}_indep")
        common.save_canvas(canvas, os.path.join(args.workdir, "10_resolution_fit" + suffix
                                                + (f"_{fixed_from_340}" if fixed_from_340 else "")
                                                + ("_allpoints" if variant == "allpoints" else "") + ".png"))

      if recipe == "uniforme" and len(fitted) >= 2:
            rows_by_resistance = {resistance: points_of(rows, resistance, fit_only=fit_only) for resistance in fitted}
            combined = fit_common(rows_by_resistance)
            if combined is None:
                continue
            values, errors = combined["result"]["values"], combined["result"]["errors"]
            canvas = ROOT.TCanvas(f"common_{variant}", "", 640 * len(fitted), 560)
            canvas.Divide(len(fitted), 1)
            print(f"  [{variant}] common fit: S {values[0]:.3f} +- {errors[0]:.3f} %, C {values[1]:.3f} +- "
                  f"{errors[1]:.3f} %, chi2/ndf {combined['result']['chi2']:.1f}/{combined['ndf']}")
            for column, resistance in enumerate(sorted(fitted)):
                index = RESISTANCES.index(resistance)
                pad = canvas.cd(column + 1)
                fit_rows = rows_by_resistance[resistance]
                draw_points(pad, fit_rows, fit_rows, resistance, selection)
                curve = common.resolution_function(common.unique_name("reso_common_curve"),
                                                   0.9 * fit_rows[0]["energy_true"], 1.05 * fit_rows[-1]["energy_true"])
                curve.SetParameters(values[2 + index], values[0], values[1])
                curve.SetLineColor(ROOT.kViolet)
                curve.SetLineStyle(2)
                curve.SetLineWidth(3)
                common.keep(curve).Draw("same")
                common.keep(common.text_box(
                    [f"N  {1000 * values[2 + index]:6.0f} #pm {1000 * errors[2 + index]:.0f} MeV",
                     f"S  {values[0]:6.3f} #pm {errors[0]:.3f} % (common)",
                     f"C  {values[1]:6.3f} #pm {errors[1]:.3f} % (common)",
                     f"#chi^{{2}}/ndf  {combined['result']['chi2']:.1f} / {combined['ndf']}"],
                    0.45, 0.66, 0.89, 0.88, 0.032, ROOT.kViolet)).Draw()
                table.append(dict(selection=selection, recipe=recipe, variant=variant, mode="common",
                                  fixed_from_340="", resistance=resistance,
                                  n_points=len(fit_rows), N_MeV=1000 * values[2 + index], err_N_MeV=1000 * errors[2 + index],
                                  N_fixed=0, S_pct=values[0], err_S=errors[0], S_fixed=0, C_pct=values[1], err_C=errors[1],
                                  C_fixed=0, chi2=combined["result"]["chi2"], ndf=combined["ndf"]))
                output.cd()
                curve.Write(f"f_{variant}_{resistance}_common")
            output.cd()
            combined["graph"].Write(f"gcomb_{variant}")
            combined["function"].Write(f"fcomb_{variant}")
            canvas.Write(f"c_{variant}_common")
            common.save_canvas(canvas, os.path.join(args.workdir, "10_resolution_common"
                                                    + ("_allpoints" if variant == "allpoints" else "") + ".png"))
    output.Close()
    common.write_csv(os.path.join(args.workdir, f"10_resolution_fits{suffix}.csv"), table, COLUMNS)


if __name__ == "__main__":
    main()
