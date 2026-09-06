#!/usr/bin/env python3
"""
Stage 11, matplotlib version -- the same final figure as s11_final_plot.py, drawn with
matplotlib from the CSVs of stages 8, 9 and 10. No ROOT needed: it runs with any python
that has numpy and matplotlib.

Writes 11_resolution[_nominal_bes]_<variant>_mpl.png, one figure per fit variant.
"""

import argparse
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESISTANCES = (340, 400, 500)
TERMS_FLOOR = 3e-4          # the terms axis does not go below this (the sync syst does)
TERMS_BY_RECIPE = {"codiceA": ("bes", "sync", "drift", "stat", "vtx_syst", "bes_syst", "syn_syst"),
                   "uniforme": ("bes", "sync", "drift", "stat", "unif_syst")}
TERM_STYLE = {"bes": ("BES", "D-.", "C1"), "sync": ("synchrotron", "^-", "C4"), "drift": ("drift", "s--", "C0"),
              "stat": ("stat", "o:", "k"), "vtx_syst": ("vtx syst", "s--", "C5"), "bes_syst": ("BES syst", "s--", "C6"),
              "syn_syst": ("sync syst", "x--", "C7"), "unif_syst": ("map syst", "P--", "C2")}
FIT_STYLE = {("indep", "S"): ("fit, S from 340 $\\Omega$", "--", "darkviolet"),
             ("indep", "SC"): ("fit, S and C from 340 $\\Omega$", ":", "green"),
             ("indep", ""): ("fit $N/E \\oplus S/\\sqrt{E} \\oplus C$", "--", "darkviolet"),
             ("common", ""): ("fit, S and C common", ":", "green")}


def read_csv(path):
    def convert(text):
        try:
            return int(text)
        except ValueError:
            try:
                return float(text)
            except ValueError:
                return text
    with open(path) as handle:
        return [{key: convert(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def resolution(energy, noise_mev, stochastic, constant):
    return np.sqrt((100 * noise_mev / 1000. / energy) ** 2 + stochastic ** 2 / energy + constant ** 2)


def parameter_text(row):
    def one(name, value, error, fixed, unit, digits):
        return f"${name}$ {value:{digits}}" + (f" {unit} (fixed)" if fixed else f" $\\pm$ {error:{digits}} {unit}")
    return "\n".join([one("N", row["N_MeV"], row["err_N_MeV"], row["N_fixed"], "MeV", "5.0f"),
                      one("S", row["S_pct"], row["err_S"], row["S_fixed"], "%", "6.3f"),
                      one("C", row["C_pct"], row["err_C"], row["C_fixed"], "%", "6.3f"),
                      f"$\\chi^2$/ndf {row['chi2']:.1f} / {row['ndf']}"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--bes", choices=("cons", "nominal"), default="cons")
    args = parser.parse_args()
    suffix = "_nominal_bes" if args.bes == "nominal" else ""
    terms_all = read_csv(os.path.join(args.workdir, "08_systematics.csv"))
    points_all = read_csv(os.path.join(args.workdir, f"09_resolution_points{suffix}.csv"))
    fits_all = [row for row in read_csv(os.path.join(args.workdir, f"10_resolution_fits{suffix}.csv"))
                if row["variant"] == "nominal"]
    selection, recipe = points_all[0]["selection"], points_all[0]["recipe"]
    resistances = [resistance for resistance in RESISTANCES if any(row["resistance"] == resistance for row in points_all)]
    bes_column = "bes_nom" if args.bes == "nominal" else "bes"
    central_label = ("$-$ BES $-$ sync" if recipe == "uniforme"
                     else f"$-$ BES ({'nominal' if args.bes == 'nominal' else 'conservative'}) $-$ sync")

    variants = []
    for row in fits_all:
        key = (row["mode"], str(row["fixed_from_340"]))
        if key not in variants:
            variants.append(key)
    for mode, fixed_from_340 in variants:
        draw_figure(args, suffix, fixed_from_340 or mode, mode, fixed_from_340, points_all, terms_all, fits_all,
                    selection, recipe, resistances, bes_column, central_label)


def draw_figure(args, suffix, name, mode, fixed_from_340, points_all, terms_all, fits_all, selection, recipe,
                resistances, bes_column, central_label):
    figure, axes = plt.subplots(2, len(resistances), figsize=(6.6 * len(resistances), 11), sharex="col",
                                gridspec_kw=dict(height_ratios=[2, 1.15]), squeeze=False)
    for column, resistance in enumerate(resistances):
        points = sorted([row for row in points_all if row["resistance"] == resistance], key=lambda row: row["energy"])
        terms = sorted([row for row in terms_all if row["resistance"] == resistance], key=lambda row: row["energy"])
        fits = [row for row in fits_all if row["resistance"] == resistance
                and (row["mode"], str(row["fixed_from_340"])) == (mode, fixed_from_340)]
        energy = np.array([row["energy_true"] for row in points], float)
        raw = np.array([row["sigma_raw"] for row in points], float)
        corrected = np.array([row["sigma_over_E"] for row in points], float)
        errors = np.array([row["err"] for row in points], float)

        axis = axes[0][column]
        axis.plot(energy, raw, "o-", ms=6, color="0.35", label="$\\sigma/\\mu$")
        axis.errorbar(energy, corrected, yerr=errors, fmt="^-", ms=7.5, color="C3", capsize=3, label=central_label)
        hand_set = np.array([str(row["window"]).endswith("-ecal_prof") for row in points])
        if hand_set.any():
            axis.plot(energy[hand_set], corrected[hand_set], "o", ms=13, mfc="none", mew=1.4, color="0.25",
                      label="no parabola: hand-set hodoscope window")
        pooled = np.array([bool(row["pooled"]) for row in points])
        if pooled.any():
            axis.plot(energy[pooled], corrected[pooled], "s", ms=14, mfc="none", mew=1.4, color="C4",
                      label="one pooled fit (no run has enough events)")
        left_out = np.array([not row["in_fit"] for row in points])
        if left_out.any():
            axis.plot(energy[left_out], corrected[left_out], "^", ms=13, mfc="none", color="0.5", label="not in the fit")
        curve_x = np.linspace(0.9 * energy.min(), 1.05 * energy.max(), 300)
        box_top = 0.62
        for fit in fits:
            label, style, colour = FIT_STYLE[(fit["mode"], str(fit["fixed_from_340"]))]
            axis.plot(curve_x, resolution(curve_x, fit["N_MeV"], fit["S_pct"], fit["C_pct"]), style, lw=2.2,
                      color=colour, label=label)
            axis.text(0.97, box_top, parameter_text(fit), transform=axis.transAxes, ha="right", va="top",
                      fontsize=9, family="monospace", bbox=dict(fc="white", ec=colour, pad=5))
            box_top -= 0.2
        axis.set_ylim(0, 1.9 * raw.max())
        axis.set_xlim(0.85 * energy.min(), 1.15 * energy.max())
        axis.set_title(f"{resistance} $\\Omega$", fontsize=12, fontweight="bold")
        axis.set_ylabel("$\\sigma/E$  [%]")
        axis.grid(alpha=.3)
        axis.legend(fontsize=8, loc="upper right")

        axis = axes[1][column]
        term_energy = np.array([row["energy_true"] for row in terms], float)
        axis.plot(term_energy, [row["sigma_raw"] for row in terms], "o-", ms=5, color="0.35", label="$\\sigma/\\mu$")
        smallest = min(row["sigma_raw"] for row in terms)
        for term in TERMS_BY_RECIPE[recipe]:
            label, style, colour = TERM_STYLE[term]
            values = np.array([row[bes_column if term == "bes" else term] for row in terms], float)
            values = np.where(np.isfinite(values) & (values > 0), values, np.nan)
            if np.isfinite(values).any():
                smallest = min(smallest, np.nanmin(values))
                axis.plot(term_energy, values, style, ms=5, color=colour, label=label)
        floor = max(0.4 * smallest, 1.6 * TERMS_FLOOR)
        single = np.array([row["drift"] <= 0 and row["n_run"] <= 1 for row in terms])
        compatible = np.array([row["drift"] <= 0 and row["n_run"] > 1 for row in terms])
        if single.any():
            axis.plot(term_energy[single], np.full(single.sum(), floor), "x", ms=6, color="C0", label="drift n/a (1 run)")
        if compatible.any():
            axis.plot(term_energy[compatible], np.full(compatible.sum(), floor), "s", ms=5, mfc="none", color="C0",
                      label="drift = 0 ($\\chi^2$/ndf $\\leq$ 1)")
        axis.set_yscale("log")
        axis.set_ylim(max(0.25 * smallest, TERMS_FLOOR), 2.5 * max(row["sigma_raw"] for row in terms))
        axis.set_xlabel("True beam energy [GeV]")
        axis.set_ylabel("size of each term  [%]")
        axis.grid(alpha=.3, which="both")
        axis.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=3, frameon=False)
    figure.suptitle(f"cut on the {selection}   $\\quad$   $A_{{tot}} > 100$ ADC", fontsize=12)
    figure.tight_layout(rect=(0, 0.02, 1, 0.98))
    path = os.path.join(args.workdir, f"11_resolution{suffix}_{name}_mpl.png")
    figure.savefig(path, dpi=150)
    print("->", path)


if __name__ == "__main__":
    main()
