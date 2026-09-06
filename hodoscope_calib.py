"""
Hodoscope calibration from the RESPONSE, so the position cut can be moved off the
ECAL centroid and onto the hodoscope.

Why not the centroid. The selection used everywhere else cuts on |pos_eta - 18| and
|pos_phi - 6|, and the centroid is built from the very amplitudes whose resolution is
being measured. Cutting on the hodoscope instead removes that circularity, but the
hodoscope reads millimetres with an arbitrary origin, so it needs an offset and a
scale first.

Which planes. There are two planes per view. In x both are healthy and their average
is taken, which is also the better position estimate. In y the two planes are equally
efficient -- both fire exactly one cluster on about 43 to 46 % of the events, y2 being
less often empty (2.4 % against 6.0 %) but correspondingly more often multi-cluster --
so efficiency does not choose between them. What chooses is the shape: profiled
against A_tot, y1 gives a clean parabola over its whole range while y2 is jagged below
zero. The plane to cut on is y1; y2 stays selectable for comparison.

How the two constants are obtained, without ever regressing against the centroid:

  offset  the crystal centre is where the response is maximum, so it is the vertex of
          the parabola fitted to <A_tot> against the hodoscope coordinate.
  scale   the same physical parabola, expressed in millimetres and in crystal units,
          has relative curvatures c_mm and c_crystal whose ratio is the square of the
          scale:  W = sqrt(c_crystal / c_mm)  millimetres per crystal. c_crystal comes
          from profili_pernorm.py.

Choosing the fit range, and knowing when there is no parabola to fit. The response is
parabolic only near the crystal centre; far out the shower leaks and the profile falls
faster, so a fit over the full range the data span puts all the lever arm on the tails
and the vertex lands wherever they pull it. No fixed range is imposed either: the beam
does not sit at the same place at every energy -- the x vertex moves from -3 mm at
20 GeV to -10 mm at 175 GeV -- and at the top energies part of the crystal falls
outside the hodoscope acceptance altogether.

So the range is not chosen, it is scanned. The profile is fitted with a parabola over
[peak - h, peak + h] for every h in SCAN_HALVES, and the answer is accepted only if it
does not depend on h:

  * at least MIN_FITS of the half-widths give a maximum inside their own fit range;
  * the vertex moves by at most VX_SPREAD_MAX over the scan;
  * W has a relative spread of at most W_SPREAD_MAX and a median inside [W_MIN, W_MAX];
  * the resulting cut window lies entirely inside the range the data span.

The vertex and W returned are the medians over the scan. Where the checks fail there
is no parabola to be found at that energy in that view, and the point is reported with
the reason rather than silently fitted: that is the hodoscope acceptance limit, and it
is drawn as such by hodo_windows.py.

Usage:
  python3 plot/hodoscope_calib.py --base . --outdir plot/hodo [--resistances 340]
"""

import argparse, os, glob, re, csv
import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import runsets

A_TOT_MIN = 100.
CORE_FRAC = 0.10                 # |A_tot/peak - 1| < CORE_FRAC for the response profile
NBINS, NMIN_BIN, NMIN_FIT = 40, 150, 3000

SCAN_HALVES = (5., 6., 7., 8., 9., 10.)   # fit half-widths scanned around the maximum [mm]
MIN_FITS = 4                              # how many of them must give a usable maximum
VX_SPREAD_MAX = 5                       # allowed excursion of the vertex over the scan [mm]
W_SPREAD_MAX = 0.7                       # allowed relative spread of W over the scan
W_MIN, W_MAX = 12., 40.                   # acceptable crystal width [mm]

FILES = {340: ("reco_340ohm", "*_merged.root"),
         400: ("reco_400ohm", "*_400_merged.root"),
         500: ("reco_500ohm", "*_500_merged.root")}


def hodo_xy(a, yplane="y1"):
    """x = mean of the two x planes, y from the requested plane; nan where unusable."""
    n1 = ak.to_numpy(a["hodo_x1_nclusters"]); n2 = ak.to_numpy(a["hodo_x2_nclusters"])
    x1 = ak.to_numpy(ak.firsts(ak.pad_none(a["hodo_x1_pos"], 1)))
    x2 = ak.to_numpy(ak.firsts(ak.pad_none(a["hodo_x2_pos"], 1)))
    ok = (n1 == 1) & (n2 == 1) & np.isfinite(x1) & np.isfinite(x2)
    x = np.where(ok, (x1 + x2) / 2., np.nan)
    ny = ak.to_numpy(a[f"hodo_{yplane}_nclusters"])
    y = ak.to_numpy(ak.firsts(ak.pad_none(a[f"hodo_{yplane}_pos"], 1)))
    return x, np.where((ny == 1) & np.isfinite(y), y, np.nan)


def response_profile(v, w, nb=NBINS, nmin=NMIN_BIN):
    """<w> against v in nb bins over the range the data span, with the error on each
    mean. No range is imposed: the profile covers the 0.5th to the 99.5th percentile.
    Returns (x, y, ey) or None."""
    m = np.isfinite(v)
    if m.sum() < NMIN_FIT:
        return None
    vv, ww = v[m], w[m]
    e = np.linspace(np.percentile(vv, 0.5), np.percentile(vv, 99.5), nb + 1)
    i = np.clip(np.digitize(vv, e) - 1, 0, nb - 1)
    xs, ys, es = [], [], []
    for k in range(nb):
        q = i == k
        if q.sum() < nmin:
            continue
        a = ww[q]
        xs.append(vv[q].mean()); ys.append(a.mean())
        es.append(a.std(ddof=1) / np.sqrt(q.sum()))
    if len(xs) < 6:
        return None
    return np.array(xs), np.array(ys), np.array(es)


def profile_peak(px, py):
    """Position of the maximum of the profile, smoothed over three bins so that a
    single fluctuating bin cannot move it."""
    sm = np.convolve(py, np.ones(3) / 3., mode="same")
    sm[0], sm[-1] = py[0], py[-1]
    return float(px[int(np.argmax(sm))])


def _fit(px, py, pe, lo, hi):
    """Weighted quadratic on the profile points inside [lo, hi]."""
    m = (px >= lo) & (px <= hi)
    n = int(m.sum())
    if n < 6:
        return None
    X = np.vstack([px[m] ** 2, px[m], np.ones(n)]).T
    iw = 1. / pe[m]
    c = np.linalg.lstsq(X * iw[:, None], py[m] * iw, rcond=None)[0]
    if c[0] >= 0:                      # no maximum: not a crystal response here
        return None
    vx = -c[1] / (2 * c[0])
    if not (lo <= vx <= hi):           # vertex extrapolated outside its own fit range
        return None
    chi2 = float((((X @ c - py[m]) / pe[m]) ** 2).sum())
    return vx, float(100 * c[0] / np.polyval(c, vx)), chi2 / max(n - 3, 1), n, np.sqrt(np.linalg.inv(X.T @ np.diag(iw**2) @ X)[1,1] / (-2*c[0])) #last is the error on vertex


def parabola_scan(px, py, pe, cc, half=0.2, halves=SCAN_HALVES):
    """Vertex and crystal width from a scan of the fit half-width.

    cc is the relative curvature of the same response in crystal units, from
    profili_pernorm.py; half is the half-window of the cut in crystal units, used only
    to check that the window fits inside the range the data span.

    Returns a dict which always carries 'ok'. When ok is False, 'why' says what went
    wrong -- that is the message to put on the plot, not a silent skip.
    """

    out = dict(ok=False, why="", x0=np.nan, W=np.nan, vx_spread=np.nan,
               w_spread=np.nan, chi2ndf=np.nan, nfit=0, peak=np.nan)
    if px is None or len(px) < 6:
        out["why"] = "no profile"
        return out
    pk = profile_peak(px, py)
    out["peak"] = pk
    if cc is None or not (cc < 0):
        out["why"] = "no curvature in crystal units"
        return out
    vs, ws, cs = [], [], []
    for h in halves:
        q = _fit(px, py, pe, pk - h, pk + h)
        if q is None:
            continue
        vs.append(q[0]); ws.append(float(np.sqrt(cc / q[1]))); cs.append(q[2])
    out["vertex_pos_error"] = q[4]
    out["nfit"] = len(vs)
    if len(vs) < MIN_FITS:
        out["why"] = f"only {len(vs)} of {len(halves)} fit ranges give a maximum"
        return out
    vs, ws = np.array(vs), np.array(ws)
    x0, W = float(np.median(vs)), float(np.median(ws))
    out.update(x0=x0, W=W, vx_spread=float(vs.max() - vs.min()),
               w_spread=float(ws.std() / W), chi2ndf=float(np.median(cs)))
    if out["vx_spread"] > VX_SPREAD_MAX:
        out["why"] = f"vertex moves by {out['vx_spread']:.1f} mm across the fit ranges"
    elif out["w_spread"] > W_SPREAD_MAX:
        out["why"] = f"W varies by {100*out['w_spread']:.0f} % across the fit ranges"
    #elif not (W_MIN <= W <= W_MAX):
    #    out["why"] = f"W = {W:.0f} mm, not a crystal"
    #elif not (px.min() <= x0 - half * W and x0 + half * W <= px.max()):
    #    out["why"] = "the window falls outside the hodoscope acceptance"
    else:
        out["ok"] = True
    return out


def crystal_curvature(plotdir, R):
    """Relative curvature in crystal units, from profili_pernorm.py."""
    out = {}
    f = os.path.join(plotdir, "profili", "profili_pernorm.csv")
    print(f)
    if not os.path.exists(f):
        print(f, " not found!!")
        return out
    for r in csv.DictReader(open(f)):
        if int(r["resistance"]) == R:
            out[(int(r["energy"]), r["coord"])] = float(r["rel_pct"])
    return out


COORD = (("x", "pos_eta"), ("y", "pos_phi"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=".")
    ap.add_argument("--outdir", default="plot/hodo")
    ap.add_argument("--plotdir", default="plot")
    ap.add_argument("--resistances", nargs="+", type=int, default=[340, 400, 500])
    ap.add_argument("--half", type=float, default=0.2)
    ap.add_argument("--yplane", choices=("y1", "y2"), default="y1")
    ap.add_argument("--exclude-runs", nargs="*", type=int, default=[])
    runsets.add_argument(ap)
    a = ap.parse_args()
    drop, only = runsets.resolve(a.runset, a.exclude_runs)
    os.makedirs(a.outdir, exist_ok=True)

    rows = []
    for R in a.resistances:
        cry = crystal_curvature(a.plotdir, R)
        d, pat = FILES[R]
        for f in sorted(glob.glob(os.path.join(a.base, d, pat)),
                        key=lambda p: int(re.match(r"(\d+)", os.path.basename(p)).group(1))):
            E = int(re.match(r"(\d+)", os.path.basename(f)).group(1))
            arr = uproot.open(f)["h4_reco"].arrays(
                ["run", "A_tot", "hodo_x1_nclusters", "hodo_x1_pos", "hodo_x2_nclusters",
                 "hodo_x2_pos", "hodo_y1_nclusters", "hodo_y1_pos",
                 "hodo_y2_nclusters", "hodo_y2_pos"], library="ak")
            run = ak.to_numpy(arr["run"]); at = ak.to_numpy(arr["A_tot"]).astype(float)
            keep = at > A_TOT_MIN
            if drop:
                keep &= ~np.isin(run, drop)
            if len(only):
                keep &= np.isin(run, only)
            if keep.sum() < NMIN_FIT:
                continue
            x, y = hodo_xy(arr, a.yplane)
            pk = float(np.median(at[keep]))
            core = keep & (np.abs(at / pk - 1) < CORE_FRAC)
            row = dict(resistance=R, energy=E, nev=int(core.sum()), peak=pk)
            msg = []
            for tag, v in (("x", x), ("y", y)):
                coord = dict(COORD)[tag]
                pr = response_profile(v[core], at[core])
                s = parabola_scan(*(pr if pr else (None, None, None)),
                                  cry.get((E, coord)), a.half)
                row.update({f"x0_{tag}": s["x0"], f"W_{tag}": s["W"],
                            f"ok_{tag}": int(s["ok"]), f"why_{tag}": s["why"],
                            f"vxspread_{tag}": s["vx_spread"],
                            f"wspread_{tag}": s["w_spread"]})
                msg.append(f"{tag}: " + (f"x0 {s['x0']:+6.2f} mm  W {s['W']:5.2f} mm"
                                         if s["ok"] else f"NO PARABOLA ({s['why']})"))
            rows.append(row)
            print(f"  {R} ohm {E:>4} GeV | " + " | ".join(msg), flush=True)

    if not rows:
        return
    cols = ("resistance,energy,nev,peak,x0_x,W_x,ok_x,vxspread_x,wspread_x,why_x,"
            "x0_y,W_y,ok_y,vxspread_y,wspread_y,why_y")
    p = os.path.join(a.outdir, "hodoscope_calib.csv")
    with open(p, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols.split(","), extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                        for k, v in r.items()})
    print("->", p)

    bad = [(r["resistance"], r["energy"], t, r[f"why_{t}"])
           for r in rows for t in "xy" if not r[f"ok_{t}"]]
    if bad:
        print("\n  no parabola at:")
        for b in bad:
            print(f"    {b[0]} ohm {b[1]:>4} GeV  {b[2]} : {b[3]}")

    fig, axs = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    for j, tag in enumerate("xy"):
        for R, mk in ((340, "o"), (400, "s"), (500, "^")):
            q = [r for r in rows if r["resistance"] == R]
            for i, key in enumerate((f"x0_{tag}", f"W_{tag}")):
                g = [r for r in q if r[f"ok_{tag}"]]
                b = [r for r in q if not r[f"ok_{tag}"] and np.isfinite(r[key])]
                if g:
                    axs[i][j].plot([r["energy"] for r in g], [r[key] for r in g], mk,
                                   ms=6, label=f"{R} $\\Omega$" if i == 0 else None)
                if b:
                    axs[i][j].plot([r["energy"] for r in b], [r[key] for r in b], mk,
                                   ms=8, mfc="none", color="C3",
                                   label="no parabola" if i == 0 and R == 340 else None)
        axs[0][j].set_ylabel(f"vertex vs hodo {tag}  [mm]")
        axs[1][j].set_ylabel(f"$W_{tag}$  [mm]")
        axs[1][j].axhline(24.2, color="0.4", lw=1, ls="--")
        axs[1][j].set_xlabel("$E_{nom}$ [GeV]")
        for i in (0, 1):
            axs[i][j].set_xscale("log"); axs[i][j].grid(alpha=.3)
        axs[0][j].legend(fontsize=8)
    fig.suptitle("hodoscope calibration from the response  "
                 "(open red = no parabola, point not used)", fontsize=12)
    fig.tight_layout()
    p = os.path.join(a.outdir, "hodoscope_calib.png")
    fig.savefig(p, dpi=150); plt.close(fig)
    print("->", p)


if __name__ == "__main__":
    main()
