"""
Response profiles against the hodoscope, with the fitted parabola and the window used.

Two figures per resistance and per view (x and y), one panel per energy: <A_tot>
normalised to its maximum against the hodoscope coordinate, the parabola accepted by
the scan in hodoscope_calib.py, its vertex, and the two vertical lines of the
selection window actually applied by resolution_hodo.py.

  hodo_windows_<view>_<R>ohm.png       axes zoomed on the parabola, one range per panel
  hodo_windows_<view>_<R>ohm_full.png  same panels on a common scale, x from -15 to
                                       +15 mm and y from 0.90 up, so the panels can be
                                       compared with each other and the tails are visible

Panels whose title is red carry no parabola: the scan found none that survives a
change of the fit range, and the reason is printed under the title. Those are the
energies where the crystal centre is at or beyond the edge of the hodoscope
acceptance, and resolution_hodo.py drops them from the parabola chain rather than
fitting them anyway. It is the picture behind the cut, and the place to look when a
point comes out odd.

Usage:
  python3 plot/hodo_windows.py --base . --outdir plot/hodo_parab \\
      --window parabola [--resistances 340]
"""

import argparse, os, glob, re
import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import runsets
from hodoscope_calib import (hodo_xy, FILES, crystal_curvature, response_profile,
                             parabola_scan, SCAN_HALVES)
from resolution_hodo import window_from_response, PLATEAU_TOL

A_TOT_MIN, CORE = 100., 0.10
LABEL = {"x": "hodo x = (x1+x2)/2  [mm]", "y": "hodo {yp}  [mm]"}
COORD = (("x", "pos_eta"), ("y", "pos_phi"))
FULL_XLIM, FULL_YLIM = (-15., 15.), (0.90, 1.005)


def draw(ax, p, cc, full):
    """One panel. p is what the reading loop collected for this energy."""
    px, py, pe, top = p["px"], p["py"], p["pe"], p["top"]
    ax.errorbar(px, py / top, yerr=pe / top, fmt="o", ms=3, lw=.8, color="C0")
    sc = p["sc"]
    # la parabola si disegna solo sull'intervallo su cui e' stata scansionata:
    # tracciarla su tutto il range darebbe l'impressione di un fit fatto anche li'
    if np.isfinite(sc["x0"]) and cc is not None and np.isfinite(sc["W"]):
        h = max(SCAN_HALVES)
        xs = np.linspace(max(sc["peak"] - h, px.min()), min(sc["peak"] + h, px.max()), 200)
        c2 = cc / (sc["W"] ** 2) / 100.
        ax.plot(xs, 1 + c2 * (xs - sc["x0"]) ** 2, "-" if sc["ok"] else "--",
                lw=1.6, color="r" if sc["ok"] else "0.6")
        ax.axvline(sc["x0"], color="r" if sc["ok"] else "0.6", lw=1, ls=":")
    if p["win"]:
        for t in p["win"]:
            ax.axvline(t, color="k", lw=1.4)
    ax.set_title(p["head"] + ("\n" + p["sub"] if p["sub"] else ""), fontsize=8.5,
                 color=p["col"], linespacing=1.5)
    ax.set_xlabel(p["xlab"], fontsize=8)
    ax.set_ylabel("$\\langle A_{tot}\\rangle$ / max", fontsize=8)
    ax.tick_params(labelsize=7); ax.grid(alpha=.3)
    if full:
        ax.set_xlim(*FULL_XLIM); ax.set_ylim(*FULL_YLIM)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=".")
    ap.add_argument("--outdir", default="plot/hodo_parab")
    ap.add_argument("--plotdir", default="plot")
    ap.add_argument("--resistances", nargs="+", type=int, default=[340, 400, 500])
    ap.add_argument("--window", choices=("plateau", "parabola"), default="parabola")
    ap.add_argument("--half", type=float, default=0.2)
    ap.add_argument("--tol", type=float, default=PLATEAU_TOL)
    ap.add_argument("--yplane", choices=("y1", "y2"), default="y1")
    ap.add_argument("--exclude-runs", nargs="*", type=int, default=[])
    runsets.add_argument(ap)
    a = ap.parse_args()
    drop, only = runsets.resolve(a.runset, a.exclude_runs)
    os.makedirs(a.outdir, exist_ok=True)

    for R in a.resistances:
        cry = crystal_curvature(a.plotdir, R)
        d, pat = FILES[R]
        files = sorted(glob.glob(os.path.join(a.base, d, pat)),
                       key=lambda p: int(re.match(r"(\d+)", os.path.basename(p)).group(1)))
        panels = {"x": [], "y": []}
        # ogni file si legge UNA volta sola e serve entrambe le viste
        for f in files:
            E = int(re.match(r"(\d+)", os.path.basename(f)).group(1))
            arr = uproot.open(f)["h4_reco"].arrays(
                ["run", "A_tot", "hodo_x1_nclusters", "hodo_x1_pos",
                 "hodo_x2_nclusters", "hodo_x2_pos", "hodo_y1_nclusters",
                 "hodo_y1_pos", "hodo_y2_nclusters", "hodo_y2_pos"], library="ak")
            run = ak.to_numpy(arr["run"]); at = ak.to_numpy(arr["A_tot"]).astype(float)
            keep = at > A_TOT_MIN
            if drop:
                keep &= ~np.isin(run, drop)
            if len(only):
                keep &= np.isin(run, only)
            if keep.sum() < 2000:
                continue
            x, y = hodo_xy(arr, a.yplane)
            core = keep & (np.abs(at / np.median(at[keep]) - 1) < CORE)
            for tag, coord in COORD:
                v = x if tag == "x" else y
                pr = response_profile(v[core], at[core])
                if pr is None:
                    continue
                px, py, pe = pr
                sc = parabola_scan(px, py, pe, cry.get((E, coord)), a.half)
                win = None
                if a.window == "parabola":
                    if sc["ok"]:
                        win = (sc["x0"] - a.half * sc["W"], sc["x0"] + a.half * sc["W"])
                else:
                    r = window_from_response(v[core], at[core], a.tol)
                    if r:
                        win = (r[0], r[1])
                sub, col = "", "k"
                if a.window == "parabola" and not sc["ok"]:
                    head, sub, col = f"{E} GeV   NO PARABOLA", sc["why"], "C3"
                elif win:
                    head = f"{E} GeV   [{win[0]:+.1f}, {win[1]:+.1f}] mm"
                    if a.window == "parabola":
                        sub = (f"$x_0$ {sc['x0']:+.2f} mm")
                else:
                    head, col = f"{E} GeV   no window", "C3"
                panels[tag].append(dict(E=E, px=px, py=py, pe=pe, top=py.max(), sc=sc,
                                        win=win, head=head, sub=sub, col=col,
                                        cc=cry.get((E, coord)),
                                        xlab=LABEL[tag].format(yp=a.yplane)))
            print(f"  {R} ohm {E:>4} GeV letto", flush=True)

        for tag, _ in COORD:
            ps = panels[tag]
            if not ps:
                continue
            nbad = sum(1 for p in ps if p["col"] == "C3")
            n = len(ps); nc = 4; nr = int(np.ceil(n / nc))
            for full in (False, True):
                fig, axs = plt.subplots(nr, nc, figsize=(4.4 * nc, 3.6 * nr),
                                        squeeze=False)
                for k, p in enumerate(ps):
                    draw(axs[k // nc][k % nc], p, p["cc"], full)
                for k in range(n, nr * nc):
                    axs[k // nc][k % nc].set_axis_off()
                head = (f"{R} $\\Omega$   $\\quad$   window: {a.window},  "
                        f"half {a.half} crystals,  fit half-widths scanned "
                        f"{SCAN_HALVES[0]:.0f}-{SCAN_HALVES[-1]:.0f} mm")
                if full:
                    head += "   $\\quad$   common scale"
                if a.window == "parabola" and nbad:
                    head += f"   $\\quad$   {nbad} energies with no parabola"
                fig.suptitle(head, fontsize=12)
                fig.tight_layout(rect=(0, 0, 1, 0.985))
                p = os.path.join(a.outdir, f"hodo_windows_{tag}_{R}ohm"
                                           f"{'_full' if full else ''}.png")
                fig.savefig(p, dpi=130); plt.close(fig)
                print("->", p)


if __name__ == "__main__":
    main()
