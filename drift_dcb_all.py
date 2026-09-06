#!/usr/bin/env python3
"""
Time drift + double-sided Crystal Ball fit per run -- ALL resistances
=====================================================================
CMS ECAL TB H4 Jun2026. Runs over every <E>[_<R>]_merged.root file found in
reco_340ohm/, reco_400ohm/, reco_500ohm/ and writes the plots to plot/<R>/.

For each (resistance, energy) it produces:
  drift_Atot_vs_evento_<E>GeV_<R>ohm.png   A_tot vs ordered event number + median profile
  dcb_fits_per_run_<E>GeV_<R>ohm.png       the individual double-CB fits, one per run
  picco_sigma_vs_run_<E>GeV_<R>ohm.png     peak, sigma and sigma/mu vs run
For each resistance:
  drift_per_run_<R>ohm.csv                 all the fit results
  sommario_drift_<R>ohm.png                run-to-run spread of the peak vs energy

Why the events have to be reordered: hadd does not mix the runs (the `run` branch is
already monotonic), but within each run the spills come in scattered order, and
within each spill the `evt` counter is rotated -- it starts from a high value and
wraps around. The time_* branches are all zero, so there is no per-event timestamp.
Chronological order is restored with a lexsort on (run, spill, evt).

Usage:
  python3 drift_dcb_all.py --base <directory with reco_*ohm/> --outdir plot \
      [--resistances 340 400 500] [--timestamps timestamps_runs.txt]
"""

import argparse
import json
import glob
import os
import re

import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import runsets
from iminuit import Minuit
from iminuit.cost import LeastSquares

# ADC->GeV scale used in fit.sh to choose the initial window
SCALE = {340: 3500 / 150., 400: 1080 / 40., 500: 3340 / 100.}

# binning identico a fit.sh:  A_tot>>h(8000, 0, 8000)
NBINS, XLO, XHI = 8000, 0., 8000.

CUT_LABEL = ("$|\\mathrm{pos\\_eta}-18| \\leq 0.2$,  "
             "$|\\mathrm{pos\\_phi}-6| \\leq 0.2$,  $A_{tot} > 0$")
BOX = ("$|\\mathrm{pos\\_eta}-18| < 0.6$,  $|\\mathrm{pos\\_phi}-6| < 0.6$")

# Noise threshold for maps and profiles:
#   base 80 ADC up to 50 GeV, 200 ADC above;
#   plus 5% of the nominal value when the latter exceeds 80 ADC.
A_TOT_BASE_LOW, A_TOT_BASE_HIGH, E_LOW = 80., 200., 50.
A_TOT_FRAC = 0.05


def a_tot_min(energy, resistance):
    base = A_TOT_BASE_LOW if energy <= E_LOW else A_TOT_BASE_HIGH
    frac = A_TOT_FRAC * SCALE[resistance] * energy
    return max(base, frac) if frac > A_TOT_BASE_LOW else base

RUN_TIME = {}          # riempito da load_timestamps()


def load_timestamps(path):
    """timestamps_runs.txt: lines of the form 'Jun 12 00:37 20769'."""
    if not path or not os.path.exists(path):
        print(f"  [!] timestamps non trovati ({path}): i plot avranno solo il numero di run")
        return
    for line in open(path):
        p = line.split()
        if len(p) >= 4 and p[-1].isdigit():
            RUN_TIME[int(p[-1])] = " ".join(p[:-1])


def position_cut(pos_eta, pos_phi):
    """taglio di fit.sh"""
    return (np.abs(pos_eta - 18) <= 0.2) & (np.abs(pos_phi - 6) <= 0.2)


# ------------------------------------------------------------ funzione DCB
def dcb_func(x, alpha_l, alpha_h, n_l, n_h, mean, sigma, N):
    """Double-sided Crystal Ball, identica a DoubleSidedCrystalballFunction in dcb.cxx."""
    t = (x - mean) / sigma
    out = np.empty_like(t)
    core = (t >= -alpha_l) & (t <= alpha_h)
    low, high = t < -alpha_l, t > alpha_h
    out[core] = np.exp(-0.5 * t[core] ** 2)
    if np.any(low):
        f2 = (n_l / alpha_l) - alpha_l - t[low]
        out[low] = np.exp(-0.5 * alpha_l ** 2) * np.power(
            np.maximum(alpha_l / n_l * f2, 1e-12), -n_l)
    if np.any(high):
        f2 = (n_h / alpha_h) - alpha_h + t[high]
        out[high] = np.exp(-0.5 * alpha_h ** 2) * np.power(
            np.maximum(alpha_h / n_h * f2, 1e-12), -n_h)
    return N * out


def hist_stats(counts, centers, lo, hi):
    """Mean and RMS of the TH1 restricted to [lo, hi] (like GetMean/GetRMS with SetRangeUser)."""
    m = (centers >= lo) & (centers <= hi)
    c, x = counts[m], centers[m]
    tot = c.sum()
    if tot <= 0:
        return np.nan, np.nan, 0.
    mean = (c * x).sum() / tot
    return mean, np.sqrt(max((c * (x - mean) ** 2).sum() / tot, 0.)), tot


def mode_window(values, energy, resistance):
    """Fallback window centred on the MODE of A_tot, for the runs whose peak does not
    sit where the fit.sh recipe puts it (e.g. 340 ohm 275 GeV, runs 20636-20639)."""
    nominal = SCALE[resistance] * energy
    v = values[(values > 0.5 * nominal) & (values < 1.3 * nominal)]
    if len(v) < 100:
        return None
    c, e = np.histogram(v, bins=150)
    mode = 0.5 * (e[c.argmax()] + e[c.argmax() + 1])
    core = v[np.abs(v - mode) < 0.08 * nominal]
    if len(core) < 50 or core.std() <= 0:
        return None
    return mode - 3 * core.std(), mode + 3 * core.std()


def fit_dcb(values, energy, resistance, n_iter=3, rebin=1, init_window=None):
    """The fit.sh recipe: window scale*E*(0.95,1.05) -> mean+-3RMS (x2) -> dcb (x3).
    With init_window the fit starts from a different window (see mode_window)."""
    nb = NBINS // rebin
    counts, edges = np.histogram(values, bins=nb, range=(XLO, XHI))
    counts = counts.astype(float)
    centers = 0.5 * (edges[:-1] + edges[1:])

    sc = SCALE[resistance]
    lo, hi = init_window if init_window else (sc * energy * 0.95, sc * energy * 1.05)
    for _ in range(2):
        mean, rms, tot = hist_stats(counts, centers, lo, hi)
        if not np.isfinite(mean) or rms <= 0:
            return None
        lo, hi = mean - 3 * rms, mean + 3 * rms
    mean, rms, tot = hist_stats(counts, centers, lo, hi)
    if tot < 50:
        return None

    sel = (centers >= lo) & (centers <= hi) & (counts > 0)   # ROOT ignora i bin vuoti
    x, y = centers[sel], counts[sel]
    if len(x) < 10:
        return None
    ey = np.sqrt(y)                                          # errori poissoniani, come TH1::Fit

    seed = dict(alpha_l=2., alpha_h=2., n_l=2., n_h=2.,
                mean=mean, sigma=rms, N=float(y.max()))
    best = None
    for _ in range(n_iter):
        m = Minuit(LeastSquares(x, y, ey, dcb_func), **seed)
        m.limits["alpha_l"] = (0.1, 10)                      # limiti identici a dcb.cxx
        m.limits["alpha_h"] = (0.1, 10)
        m.limits["n_l"] = (1, 10)
        m.limits["n_h"] = (1, 10)
        m.limits["mean"] = (lo, hi)
        m.limits["sigma"] = (0, hi - lo)
        m.limits["N"] = (0, None)
        m.migrad()
        m.hesse()
        best = m
        seed = {p: m.values[p] for p in seed}                 # re-seed from the previous fit

    # HESSE can fail ("covariance not pos. def.") when the tail parameters are
    # degenerate: with a mean+-3RMS window there are often NO bins beyond alpha_h, so
    # alpha_h/n_h do not affect the chi2 (gradient exactly 0) and the matrix is
    # singular -> all errors null. The same happens if n_l ends up on the limit of 10
    # from dcb.cxx. In those cases the tails are frozen at the fitted value and only
    # (mean, sigma, N) are refitted: peak and sigma do not change, but the errors
    # become computable (they are conditional on the shape of the tails).
    TAILS = ("alpha_l", "alpha_h", "n_l", "n_h")
    hesse_ko = (best.covariance is None or best.errors["sigma"] <= 0
                or best.errors["mean"] <= 0)
    fixed = []
    if hesse_ko:
        m = Minuit(LeastSquares(x, y, ey, dcb_func), **{p: best.values[p] for p in seed})
        for pname in TAILS:
            m.limits[pname] = best.limits[pname]
            m.fixed[pname] = True
        m.limits["mean"] = (lo, hi)
        m.limits["sigma"] = (0, hi - lo)
        m.limits["N"] = (0, None)
        m.migrad()
        m.hesse()
        if m.errors["sigma"] > 0 and m.errors["mean"] > 0:
            best, fixed = m, list(TAILS)

    mean_at_edge = (abs(best.values["mean"] - lo) < 0.02 * (hi - lo)
                    or abs(best.values["mean"] - hi) < 0.02 * (hi - lo))
    return dict(minuit=best, x=x, y=y, ey=ey, lo=lo, hi=hi, nsel=int(tot),
                chi2=best.fval, ndf=max(len(x) - 7, 1),
                peak=best.values["mean"], err_peak=best.errors["mean"],
                sigma=best.values["sigma"], err_sigma=best.errors["sigma"],
                valid=bool(best.valid) and not mean_at_edge,
                mean_at_edge=mean_at_edge, fixed="code" if fixed else "")


def _healthy(r, rb):
    if r is None:
        return False
    binw = (XHI - XLO) / (NBINS // rb)
    return (r["valid"] and r["sigma"] > 0 and r["err_sigma"] > 0
            and r["err_sigma"] / r["sigma"] < 0.25 and r["sigma"] > 2 * binw)


def fit_dcb_auto(values, energy, resistance, rebins=(1, 5, 10, 20, 40)):
    """Progressively coarser binning; if the fit.sh window misses the peak, it falls
    back to the window centred on the mode. Returns (result, rebin, window)."""
    last = (None, 0, "std")
    mw = mode_window(values, energy, resistance)
    for tag, win in (("std", None), ("mode", mw)):
        if tag == "mode" and win is None:
            continue
        for rb in rebins:
            r = fit_dcb(values, energy, resistance, rebin=rb, init_window=win)
            if r is None:
                continue
            if last[0] is None or _healthy(r, rb):
                last = (r, rb, tag)
            if _healthy(r, rb):
                return r, rb, tag
    return last



# ================================================================= sistematica
def syst_for_unit_chi2(v, e):
    """Extra error s (in quadrature) such that, fitting the points with a constant,
    chi2/ndf = 1. This is the run-to-run drift systematic.
        chi2(s) = sum (v_i - weighted_mean(s))^2 / (e_i^2 + s^2),  ndf = N-1
    If the points are already compatible (chi2/ndf <= 1) it returns s = 0."""
    v, e = np.asarray(v, float), np.asarray(e, float)
    n = len(v)
    if n < 2 or not np.all(np.isfinite(v)) or not np.all(e > 0):
        return np.nan, np.nan, np.nan

    def chi2ndf(s):
        w = 1. / (e ** 2 + s ** 2)
        m = (v * w).sum() / w.sum()
        return ((v - m) ** 2 * w).sum() / (n - 1)

    c0 = chi2ndf(0.)
    if c0 <= 1:
        w = 1. / e ** 2
        return 0., c0, (v * w).sum() / w.sum()
    lo, hi = 0., max(v.max() - v.min(), e.max())
    for _ in range(60):
        if chi2ndf(hi) <= 1:
            break
        hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if chi2ndf(mid) > 1:
            lo = mid
        else:
            hi = mid
    s = 0.5 * (lo + hi)
    w = 1. / (e ** 2 + s ** 2)
    return s, c0, (v * w).sum() / w.sum()


# ================================================== mappe 2D in (pos_eta, pos_phi)
ETA0, PHI0, HALF, NB2D = 18., 6., 0.6, 200   # griglia 200x200 di 0.006 cristalli (mappe)
NB2D_SHIFT = 12      # griglia grossa (0.1 cristalli) usata SOLO per lo shift per run:
                     # with fine bins there would be less than one event per bin
MAP_RANGE = [[ETA0 - HALF, ETA0 + HALF], [PHI0 - HALF, PHI0 + HALF]]
NMIN_RUN, NMIN_ALL = 25, 150                 # eventi minimi per bin (griglia grossa)
NMIN_MAP = 3                                 # eventi minimi per bin (griglia fine)
PAL_HALF = 0.03      # semi-ampiezza della palette attorno all'energia inclusiva


def _map(eta, phi, atot, nb=None):
    """Occupancy, <A_tot> per bin and the error on the mean per bin."""
    nb = nb or NB2D
    H, xe, ye = np.histogram2d(eta, phi, bins=nb, range=MAP_RANGE)
    S, _, _ = np.histogram2d(eta, phi, bins=nb, range=MAP_RANGE, weights=atot)
    S2, _, _ = np.histogram2d(eta, phi, bins=nb, range=MAP_RANGE, weights=atot ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(H > 0, S / np.maximum(H, 1), np.nan)
        var = np.where(H > 1, S2 / np.maximum(H, 1) - mean ** 2, np.nan)
        err = np.sqrt(np.maximum(var, 0) / np.maximum(H, 1))
    return H, mean, err, xe, ye


def centroid_figure(eta, phi, atot, run, runs, bounds, rows, energy, resistance,
                    outdir, all_eta=None, all_phi=None):
    """2D centroid, per-bin normalised <A_tot> map, and the per-run energy shift at
    fixed position (with the occupancy divided out)."""
    Hall, Mall, Eall, xe, ye = _map(eta, phi, atot)                 # griglia fine, per le mappe
    Gall, GMall, GEall, _, _ = _map(eta, phi, atot, nb=NB2D_SHIFT)   # griglia grossa, per lo shift
    cen_eta, cen_phi = eta.mean(), phi.mean()
    # centroide pesato con A_tot
    weta = (eta * atot).sum() / atot.sum()
    wphi = (phi * atot).sum() / atot.sum()

    corr = []                                  # shift per run corretto per occupancy
    for i, r in enumerate(runs):
        sl = slice(bounds[i], bounds[i + 1])
        Hr, Mr, Er, _, _ = _map(eta[sl], phi[sl], atot[sl], nb=NB2D_SHIFT)
        good = (Hr >= NMIN_RUN) & (Gall >= NMIN_ALL) & np.isfinite(Mr) & np.isfinite(GMall)
        good &= (Er > 0) & (GEall > 0)
        if good.sum() < 3:
            corr.append((np.nan, np.nan, int(good.sum())))
            continue
        ratio = Mr[good] / GMall[good]
        er = ratio * np.sqrt((Er[good] / Mr[good]) ** 2 + (GEall[good] / GMall[good]) ** 2)
        w = 1. / er ** 2
        m = (ratio * w).sum() / w.sum()
        corr.append((100 * (m - 1), 100 / np.sqrt(w.sum()), int(good.sum())))

    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1])

    ax = fig.add_subplot(gs[0, 0])
    oe = all_eta if all_eta is not None else eta
    op = all_phi if all_phi is not None else phi
    Hocc, _, _ = np.histogram2d(oe, op, bins=NB2D, range=MAP_RANGE)
    im = ax.pcolormesh(xe, ye, Hocc.T, cmap="viridis")
    fig.colorbar(im, ax=ax, label="events / bin")
    oce, ocp = oe.mean(), op.mean()
    ax.plot(oce, ocp, "rx", ms=12, mew=2.5,
            label=f"centroid ({oce:.4f}, {ocp:.4f})")
    ax.plot(weta, wphi, "w+", ms=14, mew=2.5,
            label=f"$A_{{tot}}$-weighted, above threshold ({weta:.4f}, {wphi:.4f})")
    ax.set_xlabel("pos_eta"); ax.set_ylabel("pos_phi")
    ax.set_title(f"$A_{{tot}} > 0$,   $N = {len(oe)}$", fontsize=10)
    ax.legend(fontsize=7, loc="lower left")

    ax = fig.add_subplot(gs[0, 1])
    M = np.where(Hall >= NMIN_MAP, Mall, np.nan)
    incl = atot.mean()                       # energia inclusiva
    vmin, vmax = incl * (1 - PAL_HALF), incl * (1 + PAL_HALF)
    im = ax.pcolormesh(xe, ye, M.T, cmap="plasma", vmin=vmin, vmax=vmax)
    fig.colorbar(im, ax=ax, extend="both",
                 label=f"$\\langle A_{{tot}} \\rangle$ [ADC]  "
                       f"(scale: inclusive {incl:.0f} $\\pm$ {100*PAL_HALF:.0f}%)")
    ax.plot(cen_eta, cen_phi, "kx", ms=12, mew=2.5)
    ax.set_xlabel("pos_eta"); ax.set_ylabel("pos_phi")
    ax.set_title(f"$A_{{tot}} > {a_tot_min(energy, resistance):.0f}$ ADC,   "
                 f"$N_{{bin}} \\geq {NMIN_MAP}$", fontsize=10)

    xs = np.arange(len(runs))
    ce = np.array([eta[bounds[i]:bounds[i+1]].mean() if bounds[i+1] > bounds[i] else np.nan
                   for i in range(len(runs))])
    cp = np.array([phi[bounds[i]:bounds[i+1]].mean() if bounds[i+1] > bounds[i] else np.nan
                   for i in range(len(runs))])
    we = np.array([eta[bounds[i]:bounds[i+1]].std() if bounds[i+1] - bounds[i] > 5 else np.nan
                   for i in range(len(runs))])
    wp = np.array([phi[bounds[i]:bounds[i+1]].std() if bounds[i+1] - bounds[i] > 5 else np.nan
                   for i in range(len(runs))])
    ne = np.array([max(bounds[i+1] - bounds[i], 1) for i in range(len(runs))])
    ee, ep = we / np.sqrt(ne), wp / np.sqrt(ne)

    ax = fig.add_subplot(gs[1, 0])
    ax.errorbar(xs, ce - cen_eta, yerr=ee, fmt="o-", color="C0", capsize=3, label="pos_eta")
    ax.errorbar(xs, cp - cen_phi, yerr=ep, fmt="s-", color="C1", capsize=3, label="pos_phi")
    ax.axhline(0, color="grey", lw=.8)
    ax.set_xticks(xs); ax.set_xticklabels([str(r) for r in runs], fontsize=7, rotation=90)
    ax.set_ylabel("run centroid - global centroid\n[crystal units]")
    ax.set_xlabel("run (chronological order)")

    ax.legend(fontsize=8); ax.grid(alpha=.3)

    ax = fig.add_subplot(gs[1, 1])
    pk = np.array([r["peak"] if r["ok"] else np.nan for r in rows])
    epk = np.array([r["err_peak"] if r["ok"] else np.nan for r in rows])
    p0 = np.nanmean(pk)
    raw = 100 * (pk / p0 - 1)
    eraw = 100 * epk / p0
    cshift = np.array([c[0] for c in corr], dtype=float)
    ecshift = np.array([c[1] for c in corr], dtype=float)
    # both series referred to their own mean, so that the two effects can be
    # compared directly, with no normalisation offset
    if np.isfinite(cshift).any():
        cshift = cshift - np.nanmean(cshift)
    ax.errorbar(xs, raw, yerr=eraw, fmt="o", color="C0", capsize=3,
                label="raw peak shift (double-CB fit)")
    ax.errorbar(xs + .15, cshift, yerr=ecshift, fmt="s", color="C3", capsize=3,
                label="shift at fixed position (occupancy divided out)")
    ax.axhline(0, color="grey", lw=.8)
    ax.set_xticks(xs); ax.set_xticklabels([str(r) for r in runs], fontsize=7, rotation=90)
    ax.set_ylabel("shift w.r.t. the run average [%]")
    ax.set_xlabel("run (chronological order)")

    ax.legend(fontsize=8); ax.grid(alpha=.3)

    fig.suptitle(f"{resistance} $\\Omega$,  {energy} GeV   $\\quad$   "
                 f"$A_{{tot}} > {a_tot_min(energy, resistance):.0f}$ ADC,  {BOX}",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"centroide2D_{energy}GeV_{resistance}ohm.png"), dpi=140)
    plt.close(fig)
    return [(int(r), c[0], c[1], c[2], float(a), float(b))
            for r, c, a, b in zip(runs, corr, we, wp)]



# ============ profilo di <A_tot> vs centroid (taglio sull'altra coordinata) ==
PROF_HALF, PROF_BIN, PROF_NMIN = 0.6, 0.0125, 20   # range, larghezza bin, eventi/bin minimi
SEL_HALF = 0.2       # finestra di selezione in posizione
FIT_HALF = 0.3       # intervallo su cui si fitta la parabola


def _parab(x, p0, p1, p2):
    """Parabola in the same parametrisation as the ROOT fit:
    p0 = vertex position, p1 = value at the maximum, p2 = curvature."""
    return p1 + p2 * (x - p0) ** 2


def profilo_centroide(eta, phi, atot, energy, resistance, outdir, ywin=None, amin=None):
    """<A_tot> vs (pos_eta - 18) with |pos_phi - 6| < 0.2, and vice versa.
    Profile with points and error on the mean, zoomed around the maximum, with a
    parabola overlaid (p0 = vertex, p1 = maximum, p2 = curvature).
    No constant fit."""
    nb = int(2 * PROF_HALF / PROF_BIN)
    fig, axs = plt.subplots(1, 2, figsize=(16, 6.4))
    res = []
    for ax, (x, other, xl, ol) in zip(axs, (
            (eta - ETA0, phi - PHI0, "pos_eta - 18", "pos_phi - 6"),
            (phi - PHI0, eta - ETA0, "pos_phi - 6", "pos_eta - 18"))):
        m = np.abs(other) < SEL_HALF
        if ywin:
            m = m & (atot > ywin[0]) & (atot < ywin[1])
        xx, yy = x[m], atot[m]
        H, edges = np.histogram(xx, bins=nb, range=(-PROF_HALF, PROF_HALF))
        S, _ = np.histogram(xx, bins=nb, range=(-PROF_HALF, PROF_HALF), weights=yy)
        S2, _ = np.histogram(xx, bins=nb, range=(-PROF_HALF, PROF_HALF), weights=yy ** 2)
        c = 0.5 * (edges[:-1] + edges[1:])
        ok = H >= PROF_NMIN
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = S / np.maximum(H, 1)
            var = S2 / np.maximum(H, 1) - mean ** 2
            err = np.sqrt(np.maximum(var, 0) / np.maximum(H, 1))
        ax.errorbar(c[ok], mean[ok], yerr=err[ok], fmt="o", ms=3.4, lw=.9,
                    color="C0", capsize=2, zorder=2)

        fitm = ok & (err > 0) & (np.abs(c) <= FIT_HALF)
        if fitm.sum() >= 5:
            mi = Minuit(LeastSquares(c[fitm], mean[fitm], err[fitm], _parab),
                        p0=0., p1=float(mean[fitm].max()), p2=-100.)
            mi.migrad(); mi.hesse()
            ndf = int(fitm.sum()) - 3
            xs = np.linspace(c[fitm].min(), c[fitm].max(), 400)
            ax.plot(xs, _parab(xs, *mi.values), "r-", lw=1.8, zorder=3)
            ax.text(.5, .03,
                    f"Entries {int(H[ok].sum())}    fit over |x| < {FIT_HALF}    "
                    f"$\\chi^2$/ndf {mi.fval:.1f} / {ndf}\n"
                    f"p0 (vertex) {mi.values['p0']:+.4f} $\\pm$ {mi.errors['p0']:.4f}    "
                    f"p1 (max) {mi.values['p1']:.1f} $\\pm$ {mi.errors['p1']:.1f}    "
                    f"p2 (curv.) {mi.values['p2']:.1f} $\\pm$ {mi.errors['p2']:.1f}",
                    transform=ax.transAxes, ha="center", va="bottom", fontsize=8.5,
                    bbox=dict(fc="w", ec="0.6", alpha=.92))
            res.append(dict(energy=energy, coord=xl.split()[0], nbin=int(fitm.sum()),
                            p0=mi.values["p0"], ep0=mi.errors["p0"],
                            p1=mi.values["p1"], ep1=mi.errors["p1"],
                            p2=mi.values["p2"], ep2=mi.errors["p2"],
                            chi2=mi.fval, ndf=ndf))
        # zoom attorno al massimo
        if ok.sum() > 3:
            lo_, hi_ = mean[ok].min(), mean[ok].max()
            pad = 0.12 * (hi_ - lo_) + 2 * np.median(err[ok])
            ax.set_ylim(lo_ - pad, hi_ + pad)
        for v in (-SEL_HALF, SEL_HALF):
            ax.axvline(v, color="k", lw=1)
        ax.set_xlim(-PROF_HALF, PROF_HALF)
        ax.set_xlabel(f"{xl}  [crystal units]")
        ax.set_ylabel("$\\langle A_{tot} \\rangle$ [ADC]")
        ax.set_title(f"$|\\mathrm{{{ol.replace(chr(95), chr(92)+chr(95))}}}| < {SEL_HALF}$,   "
                     f"${ywin[0]:.0f} < A_{{tot}} < {ywin[1]:.0f}$ ADC", fontsize=9.5)
        ax.grid(alpha=.3)
    fig.suptitle(f"{resistance} $\\Omega$,  {energy} GeV", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"profilo_centroide_{energy}GeV_{resistance}ohm.png"), dpi=150)
    plt.close(fig)
    return res


def scatter_centroide(eta, phi, atot, energy, resistance, outdir, ywin=None, amin=None):
    """A_tot vs (pos_eta - 18) with |pos_phi - 6| < 0.2, and vice versa.
    ALL events, no averaging and no fit: a 2D histogram with a logarithmic colour
    scale. The horizontal lines mark the effective threshold and the window in y;
    the vertical ones the +-0.2 in position."""
    nb = int(2 * PROF_HALF / PROF_BIN)
    nom = SCALE[resistance] * energy
    fig, axs2 = plt.subplots(2, 2, figsize=(16, 12))
    coords = ((eta - ETA0, phi - PHI0, "pos_eta - 18", "pos_phi - 6"),
              (phi - PHI0, eta - ETA0, "pos_phi - 6", "pos_eta - 18"))
    for row, zoom in enumerate((False, True)):
      for ax, (x, other, xl, ol) in zip(axs2[row], coords):
        m = np.abs(other) < SEL_HALF
        yr = ([ywin[0] - 3 * (ywin[1] - ywin[0]) / 20, ywin[1] + 3 * (ywin[1] - ywin[0]) / 20]
              if (zoom and ywin) else [0, 1.35 * nom])
        h = ax.hist2d(x[m], atot[m], bins=[nb, 200], range=[[-PROF_HALF, PROF_HALF], yr],
                      cmap="viridis", norm=matplotlib.colors.LogNorm(vmin=1))
        fig.colorbar(h[3], ax=ax, label="events / bin")
        if not zoom:
            ax.axhline(amin, color="red", lw=1.6,
                       label=f"threshold used for the maps: $A_{{tot}}$ > {amin:.0f} ADC "
                             f"({100 * amin / nom:.1f}% of nominal)")
        if ywin:
            for k, v in enumerate(ywin):
                ax.axhline(v, color="orange", lw=1.4, ls="--",
                           label=("window used in the 1D profile: peak $\\pm$ 10$\\sigma$ "
                                  f"({ywin[0]:.0f} - {ywin[1]:.0f} ADC)") if k == 0 else None)
        ax.set_xlabel(f"{xl}  [crystal units]")
        ax.set_ylabel("$A_{tot}$ [ADC]")
        ax.set_title(f"$|\\mathrm{{{ol.replace(chr(95), chr(92)+chr(95))}}}| < {SEL_HALF}$,   "
                     f"$A_{{tot}} > 0$", fontsize=9.5)
        ax.legend(fontsize=8, loc="lower center", framealpha=.9)
    fig.suptitle(f"{resistance} $\\Omega$,  {energy} GeV", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"Atot_vs_centroide_{energy}GeV_{resistance}ohm.png"), dpi=150)
    plt.close(fig)
    return []


def profilo_perrun(eta, phi, atot, runs, bounds, energy, resistance, outdir, ywin):
    """The same profile, but one curve per run, each normalised to its own constant:
    this shows whether the SHAPE of the response changes from run to run."""
    nb = int(2 * PROF_HALF / PROF_BIN)
    cmap = plt.get_cmap("viridis")
    fig, axs = plt.subplots(1, 2, figsize=(16, 6.4))
    for ax, (x, other, xl, ol) in zip(axs, (
            (eta - ETA0, phi - PHI0, "pos_eta - 18", "pos_phi - 6"),
            (phi - PHI0, eta - ETA0, "pos_phi - 6", "pos_eta - 18"))):
        for i, r in enumerate(runs):
            sl = slice(bounds[i], bounds[i + 1])
            m = (np.abs(other[sl]) < SEL_HALF) & (atot[sl] > ywin[0]) & (atot[sl] < ywin[1])
            if m.sum() < 500:
                continue
            xx, yy = x[sl][m], atot[sl][m]
            H, edges = np.histogram(xx, bins=nb, range=(-PROF_HALF, PROF_HALF))
            S, _ = np.histogram(xx, bins=nb, range=(-PROF_HALF, PROF_HALF), weights=yy)
            S2, _ = np.histogram(xx, bins=nb, range=(-PROF_HALF, PROF_HALF), weights=yy ** 2)
            c = 0.5 * (edges[:-1] + edges[1:])
            ok = H >= PROF_NMIN
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = S / np.maximum(H, 1)
                var = S2 / np.maximum(H, 1) - mean ** 2
                err = np.sqrt(np.maximum(var, 0) / np.maximum(H, 1))
            norm = mean[ok & (np.abs(c) <= SEL_HALF)].mean()
            ax.errorbar(c[ok], 100 * (mean[ok] / norm - 1), yerr=100 * err[ok] / norm,
                        fmt="o", ms=2.8, lw=.7, elinewidth=.7, capsize=0,
                        color=cmap(i / max(len(runs) - 1, 1)), label=str(r))
        for v in (-SEL_HALF, SEL_HALF):
            ax.axvline(v, color="k", lw=1)
        ax.axhline(0, color="grey", lw=.8)
        ax.set_xlim(-PROF_HALF, PROF_HALF)
        ax.set_xlabel(f"{xl}  [crystal units]")
        ax.set_ylabel("$\\langle A_{tot} \\rangle$ / run constant - 1  [%]")
        ax.set_title(f"$|\\mathrm{{{ol.replace(chr(95), chr(92)+chr(95))}}}| < {SEL_HALF}$,   "
                     f"${ywin[0]:.0f} < A_{{tot}} < {ywin[1]:.0f}$ ADC", fontsize=9.5)
        ax.legend(fontsize=6.5, ncol=2, title="run", title_fontsize=7)
        ax.grid(alpha=.3)
    fig.suptitle(f"{resistance} $\\Omega$,  {energy} GeV", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir,
                f"profilo_centroide_perrun_{energy}GeV_{resistance}ohm.png"), dpi=140)
    plt.close(fig)


def mappe_perrun(eta, phi, atot, runs, bounds, energy, resistance, outdir, all_eta=None,
                 all_phi=None, all_run=None, all_bounds=None):
    """For each run, the same two maps as in the centroid figure:
    occupancy (all events) and <A_tot> per bin (above threshold)."""
    ncol = min(4, len(runs))
    nrow = int(np.ceil(len(runs) / ncol))

    # --- occupancy per run, over ALL events
    if all_eta is not None and all_bounds is not None:
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.1 * ncol, 3.5 * nrow), squeeze=False)
        for ax, i, r in zip(axes.ravel(), range(len(runs)), runs):
            sl = slice(all_bounds[i], all_bounds[i + 1])
            H, xe, ye = np.histogram2d(all_eta[sl], all_phi[sl], bins=NB2D, range=MAP_RANGE)
            im = ax.pcolormesh(xe, ye, H.T, cmap="viridis")
            fig.colorbar(im, ax=ax, label="events / bin")
            ce = cp = we_ = wp_ = np.nan
            if H.sum() > 0:
                ce, cp = all_eta[sl].mean(), all_phi[sl].mean()
                ax.plot(ce, cp, "rx", ms=11, mew=2.4,
                        label=f"centroid ({ce:.4f}, {cp:.4f})")
                sl2 = slice(bounds[i], bounds[i + 1])
                if bounds[i + 1] > bounds[i] and atot[sl2].sum() > 0:
                    we_ = (eta[sl2] * atot[sl2]).sum() / atot[sl2].sum()
                    wp_ = (phi[sl2] * atot[sl2]).sum() / atot[sl2].sum()
                    ax.plot(we_, wp_, "w+", ms=13, mew=2.4,
                            label=f"pesato $A_{{tot}}$ ({we_:.4f}, {wp_:.4f})")
                ax.legend(fontsize=6.5, loc="lower left", framealpha=.85)
            ax.set_title(f"run {r},  $N = {int(H.sum())}$", fontsize=9)
            ax.set_xlabel("pos_eta", fontsize=8); ax.set_ylabel("pos_phi", fontsize=8)
            ax.tick_params(labelsize=7)
        for ax in axes.ravel()[len(runs):]:
            ax.set_axis_off()
        fig.suptitle(f"{resistance} $\\Omega$,  {energy} GeV   $\\quad$   "
                     f"$A_{{tot}} > 0$,  {BOX}", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir,
                    f"mappe2D_perrun_occupancy_{energy}GeV_{resistance}ohm.png"), dpi=130)
        plt.close(fig)

    # --- <A_tot> per bin, per run, with the same colour scale for all
    maps = []
    for i in range(len(runs)):
        sl = slice(bounds[i], bounds[i + 1])
        Hr, Mr, Er, xe, ye = _map(eta[sl], phi[sl], atot[sl])
        maps.append(np.where(Hr >= NMIN_MAP, Mr, np.nan))
    incl = atot.mean()                       # energia inclusiva
    vmin, vmax = incl * (1 - PAL_HALF), incl * (1 + PAL_HALF)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.1 * ncol, 3.5 * nrow), squeeze=False)
    for ax, i, r in zip(axes.ravel(), range(len(runs)), runs):
        im = ax.pcolormesh(xe, ye, maps[i].T, cmap="plasma", vmin=vmin, vmax=vmax)
        fig.colorbar(im, ax=ax, extend="both", label="$\\langle A_{tot} \\rangle$ [ADC]")
        ax.set_title(f"run {r},  $N = {bounds[i+1]-bounds[i]}$", fontsize=9)
        ax.set_xlabel("pos_eta", fontsize=8); ax.set_ylabel("pos_phi", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[len(runs):]:
        ax.set_axis_off()
    fig.suptitle(f"$\\langle A_{{tot}}\\rangle$ per bin, one panel per run -- "
                 f"{resistance} $\\Omega$, {energy} GeV\n"
                 f"CUTS: $A_{{tot}}$ > {a_tot_min(energy, resistance):.0f} ADC, "
                 f"|pos_eta-18| < {HALF}, |pos_phi-6| < {HALF};  same colour scale for all runs, "
                 f"narrow around the inclusive energy ({incl:.0f} $\\pm$ {100*PAL_HALF:.0f}%)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir,
                f"mappe2D_perrun_energia_{energy}GeV_{resistance}ohm.png"), dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------- analisi
def analyse(path, energy, resistance, outdir, chunk_profile, drop_runs=(), only_runs=()):

    cry_window = 3

    tree = uproot.open(path)["h4_reco"]
    arr = tree.arrays(
        ["run", "spill", "evt", "A", "pos_eta", "pos_phi",
         "sel_ieta", "sel_iphi"],
        library="ak"
    )

    arr["sel_ieta"] = ak.values_astype(arr["sel_ieta"], np.int64)
    arr["sel_iphi"] = ak.values_astype(arr["sel_iphi"], np.int64)

    if drop_runs or only_runs:
        keep_run = ~np.isin(ak.to_numpy(arr["run"]), list(drop_runs))

        if only_runs:
            keep_run &= np.isin(ak.to_numpy(arr["run"]), list(only_runs))

        if keep_run.sum() < len(keep_run):
            print(
                f"       scartati {len(keep_run)-keep_run.sum()} eventi dei run "
                f"{sorted(set(ak.to_numpy(arr['run'])[~keep_run]))}",
                flush=True
            )

        arr = arr[keep_run]

    # chronological order over ALL events (needed for the 2D maps)
    # chronological order over ALL events (needed for the 2D maps)
    o_all = np.lexsort((
        ak.to_numpy(arr["evt"]),
        ak.to_numpy(arr["spill"]),
        ak.to_numpy(arr["run"]),
    ))

    amin = a_tot_min(energy, resistance)

    # cry_window = 3 -> ieta = 17,18,19 and iphi = 5,6,7

    arr["A_total_chosen_matrix"] = ak.sum(
        arr["A"] *
        (abs(arr["sel_ieta"] - 18) < int((cry_window + 1) / 2)) *
        (abs(arr["sel_iphi"] - 6) < int((cry_window + 1) / 2)),
        axis=1
    )

    print(arr["A_total_chosen_matrix"] )

    sig = (
        (arr["A_total_chosen_matrix"][o_all] > amin)
        & (abs(arr["pos_eta"][o_all] - ETA0) < HALF)
        & (abs(arr["pos_phi"][o_all] - PHI0) < HALF)
    )

    m_run = arr["run"][o_all][sig].to_numpy()
    m_eta = arr["pos_eta"][o_all][sig].to_numpy()
    m_phi = arr["pos_phi"][o_all][sig].to_numpy()
    m_atot = arr["A_total_chosen_matrix"][o_all][sig].to_numpy()

    # same events but WITHOUT the A_tot threshold: for the occupancy and the 1D plots
    inbox = (
        (abs(arr["pos_eta"][o_all] - ETA0) < HALF)
        & (abs(arr["pos_phi"][o_all] - PHI0) < HALF)
    )

    b_eta = arr["pos_eta"][o_all][inbox].to_numpy()
    b_phi = arr["pos_phi"][o_all][inbox].to_numpy()
    b_atot = arr["A_total_chosen_matrix"][o_all][inbox].to_numpy()
    b_run = arr["run"][o_all][inbox].to_numpy()

    keep = position_cut(arr["pos_eta"].to_numpy(), arr["pos_phi"].to_numpy())

    run = arr["run"][keep].to_numpy()
    spill = arr["spill"][keep].to_numpy()
    evt = arr["evt"][keep].to_numpy()
    atot = arr["A_total_chosen_matrix"][keep].to_numpy()


    if len(atot) < 200:
        print(f"  [!] {energy} GeV: solo {len(atot)} eventi dopo il taglio, salto")
        return [], None

    # ordine cronologico vero: run -> spill -> evt (lexsort: ultima chiave = primaria)
    order = np.lexsort((evt, spill, run))
    run, atot = run[order], atot[order]
    idx = np.arange(len(atot))
    runs = np.unique(run)
    bounds = [np.searchsorted(run, r) for r in runs] + [len(run)]

    ref, _, ref_win = fit_dcb_auto(atot, energy, resistance)   # riferimento globale
    sc = SCALE[resistance]
    p0 = ref["peak"] if ref else sc * energy
    s0 = ref["sigma"] if ref else 0.01 * p0
    # -------------------------------------------------------- fit per run
    rows, fits = [], []
    for i, r in enumerate(runs):
        vals = atot[bounds[i]:bounds[i + 1]]
        res, rb, win = fit_dcb_auto(vals, energy, resistance)
        flow = float((vals < amin).mean())
        if res is None:
            rows.append(dict(run=int(r), nev=len(vals), ok=0, peak=np.nan, err_peak=np.nan,
                             sigma=np.nan, err_sigma=np.nan, chi2ndf=np.nan, rebin=0,
                             win="-", fixed="", frac_low=flow,
                             shift_occ=np.nan, err_shift_occ=np.nan, nbin_occ=0,
                             rms_eta=np.nan, rms_phi=np.nan))
            fits.append(None)
            continue
        ok = int(_healthy(res, rb))
        rows.append(dict(run=int(r), nev=len(vals), ok=ok,
                         peak=res["peak"], err_peak=res["err_peak"],
                         sigma=res["sigma"], err_sigma=res["err_sigma"],
                         chi2ndf=res["chi2"] / res["ndf"], rebin=rb,
                         win=win, fixed=res["fixed"].replace(",", "+") or "-",
                         frac_low=flow, shift_occ=np.nan, err_shift_occ=np.nan,
                         nbin_occ=0, rms_eta=np.nan, rms_phi=np.nan))
        fits.append(res if ok else None)

    # ----------------------------------------------------- pannello fit
    panel_lab = [str(r) for r in runs] + ["ALL RUNS"]
    panel_res = list(fits) + [ref]
    ncol = min(4, len(panel_lab))
    nrow = int(np.ceil(len(panel_lab) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.1 * nrow), squeeze=False)
    for ax, r, res in zip(axes.ravel(), panel_lab, panel_res):
        if res is None:
            ax.text(.5, .5, f"{r}\nfit failed", ha="center", va="center")
            ax.set_axis_off()
            continue
        ax.errorbar(res["x"], res["y"], yerr=res["ey"], fmt=".", ms=2, lw=.6, color="k")
        xs = np.linspace(res["lo"], res["hi"], 600)
        ax.plot(xs, dcb_func(xs, *[res["minuit"].values[p] for p in
                                   ("alpha_l", "alpha_h", "n_l", "n_h", "mean", "sigma", "N")]),
                "r-", lw=1.4)
        tag = (f"run {r} -- {RUN_TIME.get(int(r), '')}" if r != "ALL RUNS"
               else f"ALL RUNS together ({sum(x['nev'] for x in rows)} ev)")
        ax.set_title(f"{tag}\n"
                     f"$\\mu$={res['peak']:.1f}  $\\sigma$={res['sigma']:.2f}  "
                     f"$\\chi^2$/ndf={res['chi2']/res['ndf']:.2f}",
                     fontsize=8, fontweight=("bold" if r == "ALL RUNS" else "normal"))
        ax.set_xlabel("$A_{tot}$ [ADC]", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[len(panel_lab):]:
        ax.set_axis_off()
    fig.suptitle(f"{resistance} $\\Omega$,  {energy} GeV   $\\quad$   {CUT_LABEL}",
                 y=1.0, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"dcb_fits_per_run_{energy}GeV_{resistance}ohm.png"), dpi=140)
    plt.close(fig)

    # ------------------------------------------------------------ drift
    pk_ok = [r["peak"] for r in rows if r["ok"]]
    sg_ok = [r["sigma"] for r in rows if r["ok"]]
    if pk_ok:
        smax = max(sg_ok + [s0])
        ylo = min(pk_ok + [p0]) - 7 * smax
        yhi = max(pk_ok + [p0]) + 7 * smax
    else:
        ylo, yhi = p0 - 7 * s0, p0 + 7 * s0

    nch = max(len(idx) // chunk_profile, 1)
    cx, cy, ce = [], [], []
    for i in range(nch):
        w = atot[i * chunk_profile:(i + 1) * chunk_profile]
        w = w[(w > ylo) & (w < yhi)]
        if len(w) > 20:
            cx.append(i * chunk_profile + chunk_profile / 2)
            cy.append(np.median(w))
            ce.append(1.253 * np.std(w) / np.sqrt(len(w)))   # errore sulla mediana
    cx, cy, ce = np.array(cx), np.array(cy), np.array(ce)

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2, 1]))
    ax.hist2d(idx, atot, bins=[400, 160], range=[[0, len(idx)], [ylo, yhi]],
              cmap="viridis", cmin=1)
    if len(cx):
        ax.plot(cx, cy, "-", color="red", lw=1.5, label=f"median over {chunk_profile} events")
        ax.legend(loc="lower left", fontsize=9, framealpha=.85)
    ax.axhline(p0, color="w", ls=":", lw=1)
    ax.set_ylabel("$A_{tot}$ [ADC]")
    ax.set_ylim(ylo, yhi)
    ax.set_title(f"{resistance} $\\Omega$,  {energy} GeV,  $N = {len(idx)}$\n{CUT_LABEL}",
                 fontsize=11)

    ax2.axhline(0, color="grey", lw=1)
    if len(cx):
        ax2.errorbar(cx, 100 * (cy - p0) / p0, yerr=100 * ce / p0,
                     fmt="o", ms=3, lw=.8, color="C3", capsize=2)
    ax2.set_ylabel("(median - global peak) / peak  [%]")
    ax2.set_xlabel("event # (ordered by run $\\to$ spill $\\to$ evt)")
    ax2.grid(alpha=.3)
    for a in (ax, ax2):
        for b in bounds[1:-1]:
            a.axvline(b, color="k" if a is ax2 else "w", ls="--", lw=1, alpha=.8)
    tr = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for j, (b, r) in enumerate(zip(bounds[:-1], runs)):
        ax.text(b + 0.004 * len(idx), 0.985 - 0.11 * (j % 2),
                f"{r}\n{RUN_TIME.get(int(r), '')}", transform=tr, fontsize=7,
                va="top", ha="left", color="k",
                bbox=dict(fc="w", ec="0.6", alpha=.85, pad=1.2))
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"drift_Atot_vs_evento_{energy}GeV_{resistance}ohm.png"), dpi=150)
    plt.close(fig)


    # ---------------------------------------- centroide 2D + mappa energia
    prof = []
    if len(m_atot) > 500:
        yw = ((ref["peak"] - 10 * ref["sigma"], ref["peak"] + 10 * ref["sigma"])
              if ref else (0.8 * sc * energy, 1.2 * sc * energy))
        scatter_centroide(b_eta, b_phi, b_atot, energy, resistance, outdir,
                          ywin=yw, amin=amin)
        prof = profilo_centroide(m_eta, m_phi, m_atot, energy, resistance, outdir,
                                 ywin=yw, amin=amin)
    corr = []
    if len(m_atot) > 500:
        b2 = ([int(np.searchsorted(m_run, r, "left")) for r in runs]
              + [int(np.searchsorted(m_run, runs[-1], "right"))])
        corr = centroid_figure(m_eta, m_phi, m_atot, m_run, runs, b2, rows,
                               energy, resistance, outdir,
                               all_eta=b_eta, all_phi=b_phi)
        profilo_perrun(m_eta, m_phi, m_atot, runs, b2, energy, resistance, outdir, yw)
        ball = ([int(np.searchsorted(b_run, r, "left")) for r in runs]
                + [int(np.searchsorted(b_run, runs[-1], "right"))])
        mappe_perrun(m_eta, m_phi, m_atot, runs, b2, energy, resistance, outdir,
                     all_eta=b_eta, all_phi=b_phi, all_run=b_run, all_bounds=ball)
    cmap_shift = {c[0]: c[1:] for c in corr}
    for r in rows:
        c = cmap_shift.get(r["run"], (np.nan, np.nan, 0, np.nan, np.nan))
        (r["shift_occ"], r["err_shift_occ"], r["nbin_occ"],
         r["rms_eta"], r["rms_phi"]) = c

    # ------------------------------------------------ picco / sigma vs run
    xs = np.arange(len(rows))
    lbl = [f"{r['run']}\n{RUN_TIME.get(r['run'], '')}\n({r['nev']} ev)" for r in rows]
    okm = np.array([bool(r["ok"]) for r in rows])
    peak = np.array([r["peak"] for r in rows])
    epeak = np.array([r["err_peak"] for r in rows])
    sig = np.array([r["sigma"] for r in rows])
    esig = np.array([r["err_sigma"] for r in rows])

    fig, axs = plt.subplots(3, 1, figsize=(max(8, 1.5 * len(rows)), 10), sharex=True)
    axs[0].errorbar(xs[okm], peak[okm], yerr=epeak[okm], fmt="o", color="C0", capsize=3)
    axs[1].errorbar(xs[okm], sig[okm], yerr=esig[okm], fmt="s", color="C1", capsize=3)
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = sig / peak
        erel = rel * np.sqrt((esig / sig) ** 2 + (epeak / peak) ** 2)
    axs[2].errorbar(xs[okm], 100 * rel[okm], yerr=100 * erel[okm], fmt="^", color="C2", capsize=3)
    if ref:
        axs[0].axhline(ref["peak"], color="grey", ls="--", label=f"global fit: {ref['peak']:.1f}")
        axs[1].axhline(ref["sigma"], color="grey", ls="--", label=f"global fit: {ref['sigma']:.2f}")
        axs[2].axhline(100 * ref["sigma"] / ref["peak"], color="grey", ls="--",
                       label=f"global fit: {100*ref['sigma']/ref['peak']:.2f} %")
        for a in axs:
            a.legend(fontsize=8)
    axs[0].set_ylabel("peak $\\mu$ [ADC]")
    axs[1].set_ylabel("$\\sigma$ [ADC]")
    axs[2].set_ylabel("$\\sigma/\\mu$ [%]")
    axs[2].set_xticks(xs)
    axs[2].set_xticklabels(lbl, fontsize=7)
    axs[2].set_xlabel("run (chronological order)")
    for a in axs:
        a.grid(alpha=.3)
    axs[0].set_title(f"{resistance} $\\Omega$,  {energy} GeV\n{CUT_LABEL}", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"picco_sigma_vs_run_{energy}GeV_{resistance}ohm.png"), dpi=150)
    plt.close(fig)

    return rows, ref, prof


def summary_plot(per_energy, syst, resistance, outdir):
    """How large the run-to-run drift is, energy by energy."""
    es_all, relsig = [], []
    es_multi, spread, nrun = [], [], []
    es_single = []
    for e in sorted(per_energy):
        rows, ref = per_energy[e]
        if ref is None:
            continue
        es_all.append(e)
        relsig.append(100 * ref["sigma"] / ref["peak"])
        pk = np.array([r["peak"] for r in rows if r["ok"]])
        if len(pk) >= 2:
            es_multi.append(e)
            spread.append(100 * (pk.max() - pk.min()) / ref["peak"])
            nrun.append(len(pk))
        else:
            es_single.append(e)
    if not es_all:
        return
    fig, axs = plt.subplots(3, 1, figsize=(10, 10.5), sharex=True)
    if es_multi:
        axs[0].plot(es_multi, spread, "o", ms=8, color="C3")
        for x, y, n in zip(es_multi, spread, nrun):
            axs[0].annotate(f"{n} run", (x, y), fontsize=8, xytext=(0, 8),
                            textcoords="offset points", ha="center")
    for x in es_single:
        axs[0].annotate("1 run", (x, 0), fontsize=7, color="grey", rotation=90,
                        xytext=(0, 4), textcoords="offset points", ha="center", va="bottom")
    axs[0].axhline(0, color="grey", lw=.8)
    axs[0].set_ylabel("run-to-run peak spread\n(max-min)/peak  [%]")
    axs[0].grid(alpha=.3)
    axs[0].set_title(f"{resistance} $\\Omega$\n{CUT_LABEL}", fontsize=11)
    sy_e = [e for e in es_all if e in syst and np.isfinite(syst[e]["syst_reso_pct"])]
    if sy_e:
        axs[1].plot(sy_e, [syst[e]["syst_reso_pct"] for e in sy_e], "D-", color="C4",
                    label="systematic on resolution (sigma_E/E)")
        axs[1].legend(fontsize=8)
    axs[1].set_ylabel("drift systematic\n(error to add for $\\chi^2$/ndf=1) [%]")
    axs[1].grid(alpha=.3)

    axs[2].plot(es_all, relsig, "s-", color="C2", label="$\\sigma/\\mu$ global fit")
    axs[2].legend(fontsize=8)
    axs[2].set_ylabel("$\\sigma/\\mu$ [%]")
    axs[2].set_xlabel("Beam energy [GeV]")
    axs[2].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"sommario_drift_{resistance}ohm.png"), dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", required=True, help="cartella che contiene reco_<R>ohm/")
    p.add_argument("--outdir", default="plot")
    p.add_argument("--resistances", nargs="+", type=int, default=[340, 400, 500])
    p.add_argument("--timestamps", default=None)
    p.add_argument("--chunk-profile", type=int, default=500)
    p.add_argument("--energies", nargs="*", type=int, default=None,
                   help="rifa' solo queste energie; le altre vengono lette dalla cache")
    p.add_argument("--exclude-runs", nargs="*", type=int, default=[],
                   help="run da scartare, es. 20592 (energia sbagliata)")
    runsets.add_argument(p)
    p.add_argument("--only-collect", action="store_true",
                   help="non rifa' nessun fit: assembla CSV e sommario dalla cache")
    a = p.parse_args()
    drop, only = runsets.resolve(a.runset, a.exclude_runs)

    load_timestamps(a.timestamps or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                 "timestamps_runs.txt"))

    for R in a.resistances:
        d = os.path.join(a.base, f"reco_{R}ohm")
        files = sorted(glob.glob(os.path.join(d, "*_merged.root")))

        if not files:
            print(f"[{R} ohm] nessun file in {d}, salto")
            continue
        out = os.path.join(a.outdir, str(R))
        os.makedirs(out, exist_ok=True)
        print(f"[{R} ohm] {len(files)} file -> {out}")

        cache_dir = os.path.join(out, "_cache")
        os.makedirs(cache_dir, exist_ok=True)

        # The per-energy fits are the slow part. Each energy is cached as soon as it
        # is done, so the work can be split over several runs of the script
        # (--energies) and assembled at the end (--only-collect) with nothing redone.
        for f in files:
            m = re.match(r"^(\d+)_", os.path.basename(f))
            if not m:
                continue
            E = int(m.group(1))
            cf = os.path.join(cache_dir, f"{E}.json")
            if a.only_collect:
                continue
            if a.energies is not None and E not in a.energies:
                continue
            print(f"  {E:4d} GeV  {os.path.basename(f)}", flush=True)
            rows, ref, prof = analyse(f, E, R, out, a.chunk_profile, drop, only)
            slim = None if ref is None else {k: float(ref[k]) for k in
                                             ("peak", "err_peak", "sigma", "err_sigma")}
            if slim is not None:
                slim["chi2ndf"] = float(ref["chi2"] / ref["ndf"])
                slim["fixed"] = ref["fixed"]
            with open(cf, "w") as fh:
                json.dump(dict(rows=rows, ref=slim, prof=prof), fh, default=float)
            print(f"       -> cache {cf}", flush=True)

        per_energy, csv_rows, prof_rows = {}, [], []
        for f in files:
            m = re.match(r"^(\d+)_", os.path.basename(f))
            if not m:
                continue
            E = int(m.group(1))
            cf = os.path.join(cache_dir, f"{E}.json")
            if not os.path.exists(cf):
                print(f"  [!] {E} GeV: cache assente, salto")
                continue
            with open(cf) as fh:
                c = json.load(fh)
            rows, ref, prof = c["rows"], c["ref"], c["prof"]
            prof_rows += prof
            if not rows:
                continue
            per_energy[E] = (rows, ref)
            for r in rows:
                csv_rows.append((E, r["run"], r["nev"], r["ok"], r["rebin"], r["win"],
                                 r["fixed"], round(r["frac_low"], 4), r["peak"],
                                 r["err_peak"], r["sigma"], r["err_sigma"], r["chi2ndf"],
                                 r["shift_occ"], r["err_shift_occ"], r["nbin_occ"],
                                 r["rms_eta"], r["rms_phi"]))
            if ref:
                csv_rows.append((E, "ALL", sum(r["nev"] for r in rows), 1, 1, "std",
                                 ref["fixed"].replace(",", "+") or "-",
                                 round(float(np.mean([r["frac_low"] for r in rows])), 4),
                                 ref["peak"], ref["err_peak"], ref["sigma"],
                                 ref["err_sigma"], ref["chi2ndf"],
                                 float("nan"), float("nan"), 0,
                                 float("nan"), float("nan")))

        with open(os.path.join(out, f"drift_per_run_{R}ohm.csv"), "w") as fh:
            fh.write("energy,run,nev,fit_ok,rebin,window,par_fissati,frac_A_tot_basso,"
                     "peak_abs,err_peak_abs,sigma_abs,err_sigma_abs,chi2_ndf,"
                     "shift_occ_pct,err_shift_occ_pct,nbin_occ,rms_eta,rms_phi\n")
            for c in csv_rows:
                fh.write(",".join(f"{v:.4f}" if isinstance(v, float) else str(v) for v in c) + "\n")
        # ------------------------------------------------- sistematica di drift
        syst = {}
        for e in sorted(per_energy):
            rows, ref = per_energy[e]
            good = [r for r in rows if r["ok"]]
            if ref is None or len(good) < 2:
                print("not saved")
                continue
            pk = np.array([r["peak"] for r in good])
            epk = np.array([r["err_peak"] for r in good])
            sg = np.array([r["sigma"] for r in good])
            esg = np.array([r["err_sigma"] for r in good])
            ss, c0s, ms = syst_for_unit_chi2(sg/pk, sg/pk * ((epk/pk)**2 + (esg/sg)**2)**0.5 )
            syst[e] = dict(n=len(good), reso_w=ms, syst_reso=ss, chi2_reso=c0s,
                           syst_reso_pct=ss*100)
        with open(os.path.join(out, f"sistematica_drift_{R}ohm.csv"), "w") as fh:
            fh.write("energy,n_run,reso_medio,chi2_ndf_resp_senza_syst,syst_reso,"
                     "syst_reso_pct,reso_globale\n")
            for e in sorted(syst):
                d = syst[e]
                sg_glob = per_energy[e][1]["sigma"]/per_energy[e][1]["peak"]
                fh.write(f"{e},{d['n']},{d['reso_w']:.4f},{d['chi2_reso']:.4f},"
                         f"{d['syst_reso']:.4f},{d['syst_reso_pct']:.4f},{sg_glob:.4f}\n")
        with open(os.path.join(out, f"parabola_centroide_{R}ohm.csv"), "w") as fh:
            fh.write("energy,coord,n_bin,p0_vertice,err_p0,p1_massimo,err_p1,"
                     "p2_curvatura,err_p2,chi2,ndf\n")
            for d in prof_rows:
                fh.write(f"{d['energy']},{d['coord']},{d['nbin']},{d['p0']:.6f},"
                         f"{d['ep0']:.6f},{d['p1']:.4f},{d['ep1']:.4f},{d['p2']:.4f},"
                         f"{d['ep2']:.4f},{d['chi2']:.4f},{d['ndf']}\n")
        summary_plot(per_energy, syst, R, out)
        print(f"[{R} ohm] fatto")


if __name__ == "__main__":
    main()
