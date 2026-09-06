#!/usr/bin/env python3
"""
Stage 4 -- the response surface of every run, and save.

Reads 03_moments.csv and solves, run by run, the normal equations M c = b of the
weighted linear least squares for

    f(u, v) = a0 + a1 u + a2 v + a3 u^2 + a4 v^2 + a5 uv      u = pos_eta - 18, v = pos_phi - 6

on the events inside the fit window of the run, normalised to its peak (so a0 ~ 1 and
a3/a0, a4/a0 are the relative curvatures in %/crystal^2 once multiplied by 100). This is
the "run" variant of uniformita_maps.py: a run needs at least MIN_EVENTS_OWN_SURFACE
events in the window to have a surface of its own; below that ok = 0 and stage 6 falls
back on the surface of the energy.

The solve is done with TMatrixD (ROOT): 6x6 per run, nothing is fitted again.

Writes 04_surface_per_run.csv and surfaces/curvature_per_run_<R>ohm.png.
"""

import argparse
import os

import numpy as np
import ROOT

import common

MIN_EVENTS_OWN_SURFACE = 200     # NMIN_RUN of uniformita_maps.py
COLUMNS = (["resistance", "energy", "run", "n_window", "ok"] + [f"a{index}" for index in range(6)]
           + ["curv_eta_pct", "curv_phi_pct", "cross_pct"])


def moment_matrix(row):
    matrix = ROOT.TMatrixD(6, 6)
    vector = ROOT.TVectorD(6)
    for row_index in range(6):
        for column_index in range(6):
            matrix[row_index][column_index] = row[f"M_{row_index}{column_index}"]
        vector[row_index] = row[f"b_{row_index}"]
    return matrix, vector


def solve_surface(matrix, vector):
    """c = M^-1 b by LU decomposition, or None when M is singular."""
    solution = ROOT.TVectorD(vector)
    decomposition = ROOT.TDecompLU(matrix)
    if not decomposition.Solve(solution):
        return None
    coefficients = [float(solution[index]) for index in range(6)]
    if not all(np.isfinite(coefficients)):
        return None
    return coefficients


def curvature_columns(coefficients):
    constant = coefficients[0]
    return dict(curv_eta_pct=100 * coefficients[3] / constant,
                curv_phi_pct=100 * coefficients[4] / constant,
                cross_pct=100 * coefficients[5] / constant)


def draw_curvatures(rows, resistance, outdir):
    rows = [row for row in rows if row["resistance"] == resistance and row["ok"]]
    if not rows:
        return
    canvas = ROOT.TCanvas(f"curv_run_{resistance}", "", 1300, 500)
    canvas.Divide(2, 1)
    for pad_index, (key, label) in enumerate((("curv_eta_pct", "a3/a0 (eta)"), ("curv_phi_pct", "a4/a0 (phi)"))):
        pad = canvas.cd(pad_index + 1)
        pad.SetLogx()
        pad.SetGrid()
        graph = common.keep(common.make_graph([row["energy"] for row in rows], [row[key] for row in rows]))
        graph.SetMarkerStyle(common.MARKER[resistance])
        graph.SetMarkerColor(common.COLOUR[resistance])
        graph.SetTitle(f"{resistance} #Omega, one point per run;E_{{nom}} [GeV];{label} [% / crystal^{{2}}]")
        graph.Draw("AP")
        graph.GetXaxis().SetLimits(15, 320)
    common.save_canvas(canvas, os.path.join(outdir, "surfaces", f"curvature_per_run_{resistance}ohm.png"))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", required=True)
    args = parser.parse_args()
    common.style()

    rows = []
    for moment_row in common.read_csv(common.require(os.path.join(args.workdir, "03_moments.csv"))):
        row = dict(resistance=moment_row["resistance"], energy=moment_row["energy"], run=moment_row["run"],
                   n_window=moment_row["n_window"], ok=0)
        coefficients = None
        if moment_row["n_window"] >= MIN_EVENTS_OWN_SURFACE:
            coefficients = solve_surface(*moment_matrix(moment_row))
        if coefficients is not None:
            row["ok"] = 1
            row.update({f"a{index}": value for index, value in enumerate(coefficients)})
            row.update(curvature_columns(coefficients))
        else:
            row.update({f"a{index}": np.nan for index in range(6)})
            row.update(curv_eta_pct=np.nan, curv_phi_pct=np.nan, cross_pct=np.nan)
        rows.append(row)
    common.write_csv(os.path.join(args.workdir, "04_surface_per_run.csv"), rows, COLUMNS)
    for resistance in sorted({row["resistance"] for row in rows}):
        draw_curvatures(rows, resistance, args.workdir)
        without = [row["run"] for row in rows if row["resistance"] == resistance and not row["ok"]]
        if without:
            print(f"  {resistance} ohm: {len(without)} runs without their own surface "
                  f"(fewer than {MIN_EVENTS_OWN_SURFACE} events in the window): {without}")


if __name__ == "__main__":
    main()
