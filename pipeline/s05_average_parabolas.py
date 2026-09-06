#!/usr/bin/env python3
"""
Stage 5 -- the average response surfaces, and save.

Reads 03_moments.csv and sums the moments over runs: the surface of one ENERGY is the
solution of the summed normal equations of all its runs, the surface of one RESISTANCE
("mean") the same over all its energies. Because M and b are sufficient statistics of
the weighted least squares, summing them is exactly the fit on the union of the events
(uniformita_maps.py, coef_ene and coef_mean).

Writes 05_surface_mean.csv, one row per (resistance, energy) with scope = energy and one
per resistance with scope = mean and energy = 0, and surfaces/curvature_mean_<R>ohm.png.
"""

import argparse
import os
from collections import defaultdict

import ROOT

import common
from s04_fit_parabolas import curvature_columns, moment_matrix, solve_surface

COLUMNS = (["resistance", "energy", "scope", "n_runs", "n_window"] + [f"a{index}" for index in range(6)]
           + ["curv_eta_pct", "curv_phi_pct", "cross_pct"])


def summed(rows):
    matrix, vector = ROOT.TMatrixD(6, 6), ROOT.TVectorD(6)
    for row in rows:
        row_matrix, row_vector = moment_matrix(row)
        matrix += row_matrix
        vector += row_vector
    return matrix, vector


def surface_row(rows, resistance, energy, scope):
    coefficients = solve_surface(*summed(rows))
    row = dict(resistance=resistance, energy=energy, scope=scope, n_runs=len(rows),
               n_window=sum(item["n_window"] for item in rows))
    if coefficients is None:
        print(f"  {resistance} ohm {energy} GeV ({scope}): singular surface")
        return None
    row.update({f"a{index}": value for index, value in enumerate(coefficients)})
    row.update(curvature_columns(coefficients))
    return row


def draw(rows, resistance, outdir):
    energy_rows = sorted([row for row in rows if row["resistance"] == resistance and row["scope"] == "energy"],
                         key=lambda row: row["energy"])
    mean_rows = [row for row in rows if row["resistance"] == resistance and row["scope"] == "mean"]
    if not energy_rows:
        return
    canvas = ROOT.TCanvas(f"curv_mean_{resistance}", "", 1300, 500)
    canvas.Divide(2, 1)
    for pad_index, (key, label) in enumerate((("curv_eta_pct", "a3/a0 (eta)"), ("curv_phi_pct", "a4/a0 (phi)"))):
        pad = canvas.cd(pad_index + 1)
        pad.SetLogx()
        pad.SetGrid()
        graph = common.keep(common.make_graph([row["energy"] for row in energy_rows],
                                              [row[key] for row in energy_rows]))
        graph.SetMarkerStyle(common.MARKER[resistance])
        graph.SetMarkerColor(common.COLOUR[resistance])
        graph.SetLineColor(common.COLOUR[resistance])
        graph.SetTitle(f"{resistance} #Omega, surface per energy;E_{{nom}} [GeV];{label} [% / crystal^{{2}}]")
        graph.Draw("APL")
        graph.GetXaxis().SetLimits(15, 320)
        if mean_rows:
            line = common.keep(ROOT.TLine(15, mean_rows[0][key], 320, mean_rows[0][key]))
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray + 2)
            line.Draw()
            legend = common.keep(ROOT.TLegend(0.5, 0.78, 0.88, 0.88))
            legend.AddEntry(graph, "per energy", "pl")
            legend.AddEntry(line, "all energies (mean surface)", "l")
            legend.Draw()
    common.save_canvas(canvas, os.path.join(outdir, "surfaces", f"curvature_mean_{resistance}ohm.png"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()
    common.style()

    moment_rows = common.read_csv(common.require(os.path.join(args.workdir, "03_moments.csv")))
    by_energy, by_resistance = defaultdict(list), defaultdict(list)
    for row in moment_rows:
        by_energy[(row["resistance"], row["energy"])].append(row)
        by_resistance[row["resistance"]].append(row)

    rows = []
    for (resistance, energy), group in sorted(by_energy.items()):
        row = surface_row(group, resistance, energy, "energy")
        if row:
            rows.append(row)
            print(f"  {resistance} ohm {energy:>4} GeV: c_eta {row['curv_eta_pct']:+.2f}  c_phi "
                  f"{row['curv_phi_pct']:+.2f}  cross {row['cross_pct']:+.2f} %/crystal^2 ({len(group)} runs)")
    for resistance, group in sorted(by_resistance.items()):
        row = surface_row(group, resistance, 0, "mean")
        if row:
            rows.append(row)
            print(f"  {resistance} ohm mean surface: c_eta {row['curv_eta_pct']:+.2f}  c_phi "
                  f"{row['curv_phi_pct']:+.2f} %/crystal^2")
    common.write_csv(os.path.join(args.workdir, "05_surface_mean.csv"), rows, COLUMNS)
    for resistance in sorted(by_resistance):
        draw(rows, resistance, args.workdir)


if __name__ == "__main__":
    main()
