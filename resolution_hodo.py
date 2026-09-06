
"""
Energy resolution with the position cut made on the HODOSCOPE instead of the centroid.

The standard selection cuts on the ECAL centroid, |pos_eta - 18| <= 0.2 and
|pos_phi - 6| <= 0.2. The centroid is computed from the same amplitudes whose width
is being measured, so the cut is not independent of the quantity under study. Here
the same geometric cut is applied to the hodoscope, which measures the impact point
independently:

    |hodo_x - x0_x| <= HALF * W_x        |hodo_y - x0_y| <= HALF * W_y

with x0 and W from hodoscope_calib.py: x0 is the vertex of the response parabola in
millimetres, W the crystal width recovered from the ratio of the curvatures in
millimetres and in crystal units. HALF is 0.2, the same half-window as the centroid
cut.

Planes: x is the average of the two x planes, y is y1 -- the two y planes are equally
efficient, but y2's profile is jagged below zero while y1 is parabolic over its whole
range. Only events with exactly one cluster in EACH plane used are kept -- x1, x2 and
y1, all three -- which leaves about 35 % of the events. That cost is reported per
point, and it is the price of a cut that does not depend on the amplitudes.

Two ways of setting the window, --window:
  parabola  vertex of the response parabola +- half * W, with the fit range scanned
            and the answer accepted only if it does not depend on the range. Where no
            parabola survives the scan the energy is dropped and the reason printed.
  plateau   the contiguous range around the maximum where the response stays within
            --tol of it. No calibration at all, and no parabola needed.

From there the chain is the usual one: double-CB fit per run, weighted mean of the
per-run sigma/peak, then BES and synchrotron subtracted in quadrature.

Usage:
  python3 plot/resolution_hodo.py --base . --outdir plot/hodo --besdir plot/bes \\
      [--resistances 340] [--half 0.2]
"""

import argparse, os, glob, re, csv, math
import numpy as np
import uproot
import awkward as ak
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares
import runsets
from uniformita_pos import (fit_dcb, rel, wmean, wscatter, design, VERA,
                           A_TOT_MIN, TAILS)
from drift_dcb_all import syst_for_unit_chi2
from hodoscope_calib import (hodo_xy, FILES, crystal_curvature, response_profile,
                             parabola_scan, COORD)

SYNC_C = 1.92e-7
ETA0, PHI0, SEL = 18., 6., 0.182

PLATEAU_TOL = 0.005          # frazione sotto il massimo che conta ancora come piatta

BAD_PARABOLA_VTX_X = {(500, 50): -6.5}
BAD_PARABOLA_VTX_Y = {(400, 20): 11, (340, 225): 4, (500, 80): -6.8, (500, 30): 8, (500, 40): -5, (500, 50): 2, (500, 60): -6}
BAD_PARABOLA_WIDTH_Y = {(400, 20): 0.09}
BAD_PARABOLA_WIDTH_X = {}


def window_from_response(h, at, tol, lo=None, hi=None, nb=40, nmin=150):
    """Hodoscope window where the response is flat, without any calibration.

    The cut on the centroid selects the region around the crystal centre where the
    response varies little. The same region can be defined directly on the hodoscope,
    with no offset and no scale: profile <A_tot> against the hodoscope coordinate and
    keep the contiguous range around the maximum where the profile stays within tol of
    its plateau value. Nothing here uses the centroid, so the selection does not
    depend on the amplitudes whose width is being measured.

    No range is imposed in either view: the profile covers the whole range the data
    span, from the 0.5th to the 99.5th percentile, and the window is set only by where
    the response is flat. That is safe on y1, whose profile is a clean response curve
    over its whole range; on y2, which is jagged below zero, the plateau search would
    lock onto a fluctuation instead.

    Returns (lo, hi, drop_pct, n) with drop_pct the response variation across the
    window, which is what makes windows at different energies comparable.
    """
    m = np.isfinite(h)
    if lo is not None:
        m &= (h > lo)
    if hi is not None:
        m &= (h < hi)
    if m.sum() < 3000:
        return None
    # nessun limite a priori: il profilo si fa su tutto il range coperto dai dati e
    # la finestra la decide solo la piattezza della risposta
    e = np.linspace(np.percentile(h[m], 0.5), np.percentile(h[m], 99.5), nb + 1)
    i = np.clip(np.digitize(h[m], e) - 1, 0, nb - 1)
    xs, ys, ns = [], [], []
    for k in range(nb):
        q = i == k
        if q.sum() < nmin:
            continue
        xs.append(h[m][q].mean()); ys.append(at[m][q].mean()); ns.append(q.sum())
    if len(xs) < 6:
        return None
    xs, ys = np.array(xs), np.array(ys)
    # massimo su tre bin, cosi' una fluttuazione singola non sposta il centro
    sm = np.convolve(ys, np.ones(3) / 3., mode="same")
    sm[0], sm[-1] = ys[0], ys[-1]
    k0 = int(np.argmax(sm))
    thr = sm[k0] * (1 - tol)
    a = k0
    while a > 0 and sm[a - 1] >= thr:
        a -= 1
    b = k0
    while b < len(sm) - 1 and sm[b + 1] >= thr:
        b += 1
    if b <= a:
        return None
    drop = 100 * (1 - min(ys[a], ys[b]) / ys[k0])
    return float(xs[a]), float(xs[b]), float(drop), int(sum(ns[a:b + 1]))


def load_bes(besdir, R):
    b_nom, b_cons = {}, {}
    f = os.path.join(besdir, f"colls_energies_summary_{R}ohm.csv")
    print(f)
    if os.path.exists(f):
        for r in csv.DictReader(open(f)):
            b_nom[int(float(r["Energy"]))] = float(r["BES_formula"])
            b_cons[int(float(r["Energy"]))] = float(r["BES_cons"]) #TO FIX!!
    return b_nom, b_cons


def reso(x, N, S, C):
    return np.sqrt((100 * N / x) ** 2 + (S / np.sqrt(x)) ** 2 + C ** 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=".")
    ap.add_argument("--cry_width_mm", default=22)
    ap.add_argument("--outdir", default="plot/hodo")
    ap.add_argument("--besdir", default="plot/bes")
    ap.add_argument("--plotdir", default="plot")
    ap.add_argument("--resistances", nargs="+", type=int, default=[340])
    ap.add_argument("--half", type=float, default=SEL,
                    help="half-window of the centroid cut, in crystal units")
    ap.add_argument("--window", choices=("plateau", "parabola"), default="plateau",
                    help="plateau = range where the response stays within --tol of its "
                         "maximum, no calibration; parabola = vertex of the response "
                         "parabola +- half * W, with W the crystal width from the ratio "
                         "of the curvatures in mm and in crystal units")
    ap.add_argument("--tol", type=float, default=PLATEAU_TOL,
                    help="how far below the maximum the response may fall inside the "
                         "hodoscope window; 0.005 means the plateau is kept to 0.5%%")
    ap.add_argument("--yplane", choices=("y1", "y2"), default="y1",
                    help="hodoscope y plane to cut on; y1 is the one with a clean response profile over its whole range")
    ap.add_argument("--exclude-runs", nargs="*", type=int, default=[])
    ap.add_argument("--exclude", nargs="*", default=["340:275"],
                    help="R:E points dropped entirely")
    ap.add_argument("--nofit-energies", nargs="*", type=int, default=[],
                    help="energies kept in the plot but left out of the N/S/C fit")
    ap.add_argument("--tails", choices=("free", "fixed", "both"), default="both",
                    help="DCB tail parameters per run: free, or held at the values of "
                         "the pooled fit of the same energy. 'both' takes free as the "
                         "nominal and carries |free - fixed| as a systematic on the "
                         "fit model, added to the error bar")
    runsets.add_argument(ap)
    a = ap.parse_args()
    drop, only = runsets.resolve(a.runset, a.exclude_runs)
    excl = {tuple(int(v) for v in s.split(":")) for s in a.exclude}
    os.makedirs(a.outdir, exist_ok=True)

    rows = []
    for R in a.resistances:
        bes_nom, bes_cons = load_bes(a.besdir, R)

        cry = crystal_curvature(a.plotdir, R)
        d, pat = FILES[R]
        for f in sorted(glob.glob(os.path.join(a.base, d, pat)),
                        key=lambda p: int(re.match(r"(\d+)", os.path.basename(p)).group(1))):
            E = int(re.match(r"(\d+)", os.path.basename(f)).group(1))
            if (R, E) in excl:
                continue
            arr = uproot.open(f)["h4_reco"].arrays(
                ["run", "A", "A_tot", "sel_ieta", "sel_iphi", "pos_eta", "pos_phi", "hodo_x1_nclusters", "hodo_x1_pos",
                 "hodo_x2_nclusters", "hodo_x2_pos", "hodo_y1_nclusters", "hodo_y1_pos",
                 "hodo_y2_nclusters", "hodo_y2_pos"],
                library="ak")


            arr["sel_ieta"] = ak.values_astype(arr["sel_ieta"], np.int64)
            arr["sel_iphi"] = ak.values_astype(arr["sel_iphi"], np.int64)


            run = ak.to_numpy(arr["run"]);

            cry_window = 3

            a_3x3 = ak.to_numpy(ak.sum(
                arr["A"] *
                (abs(arr["sel_ieta"] - 18) < int((cry_window + 1) / 2) ) *
                (abs(arr["sel_iphi"] - 6) < int((cry_window + 1) / 2) ),
                axis=1
            ))

            at = ak.to_numpy(arr["A_tot"])

            pe = ak.to_numpy(arr["pos_eta"]); pp = ak.to_numpy(arr["pos_phi"])
            x, y = hodo_xy(arr, a.yplane)
            base = at > A_TOT_MIN
            if drop:
                base &= ~np.isin(run, drop)
            if len(only):
                base &= np.isin(run, only)
            cut_c = {"nominal": base & (np.abs(pe - ETA0) <= a.half) & (np.abs(pp - PHI0) <= a.half)}

            # finestra dove la risposta e' piatta, senza calibrazione: nessun offset,
            # nessuna scala, nessun uso del centroide
            core = base & (np.abs(at / np.median(at[base]) - 1) < 0.10)
            why, fallback = {}, []
            if a.window == "parabola":
                # offset dal vertice della parabola di risposta, scala dal rapporto
                # fra la curvatura in mm e quella in unita' di cristallo. Il range di
                # fit non e' fissato: viene scansionato, e la parabola si accetta solo
                # se vertice e W non dipendono dal range (vedi hodoscope_calib.py).
                wx = wy = None
                vertex_error = {"x": 1, "y": 1}
                for tag, v in (("x", x), ("y", y)):
                    pr = response_profile(v[core], at[core])
                    sc = parabola_scan(*(pr if pr else (None, None, None)),
                                       cry.get((E, dict(COORD)[tag])), a.half)
                    why[tag] = "" if sc["ok"] else sc["why"]
                    if not sc["ok"]:
                        continue
                    w = (sc["x0"] - a.half * a.cry_width_mm, sc["x0"] + a.half * a.cry_width_mm,
                         0., 0)
                    if tag == "x":
                        wx = w
                    else:
                        wy = w

                # dove la parabola non esiste il punto non si butta via: si ricade
                # sulla finestra di plateau, che non ha bisogno di nessun fit, e il
                # punto resta nel grafico marcato come tale. Buttarlo lo farebbe
                # sparire dal plot senza dire perche'
                if (R, E) in BAD_PARABOLA_VTX_X: #se va a pouttane rimettere (if wx is None)
                    print("fallback x", R, E)
                    if (R, E) in BAD_PARABOLA_WIDTH_X: width = BAD_PARABOLA_WIDTH_X[(R, E)]
                    else: width = a.half
                    wx = (BAD_PARABOLA_VTX_X[(R, E)] - width * a.cry_width_mm, BAD_PARABOLA_VTX_X[(R, E)] + width * a.cry_width_mm, 0, 0)
                    fallback.append("x")
                if (R, E) in BAD_PARABOLA_VTX_Y:
                    print("fallback y", R, E)
                    if (R, E) in BAD_PARABOLA_WIDTH_Y: width = BAD_PARABOLA_WIDTH_Y[(R, E)]
                    else: width = a.half
                    wy = (BAD_PARABOLA_VTX_Y[(R, E)] - width * a.cry_width_mm, BAD_PARABOLA_VTX_Y[(R, E)] + width * a.cry_width_mm, 0, 0)
                    fallback.append("y")

            print(fallback)

            if fallback:
                print(f"  {R} ohm {E:>4} GeV: niente parabola in "
                      f"{'+'.join(fallback)} ({'; '.join(why[t] for t in fallback)}), "
                      f"finestra dal plateau")

            cut_h = {
              "nominal": (base & np.isfinite(x) & np.isfinite(y) & (x >= wx[0]) & (x <= wx[1]) & (y >= wy[0]) & (y <= wy[1])),
              "x_low": (base & np.isfinite(x) & np.isfinite(y) & (x >= wx[0] - vertex_error["x"]) & (x <= wx[1] - vertex_error["x"]) & (y >= wy[0]) & (y <= wy[1])),
              "x_high": (base & np.isfinite(x) & np.isfinite(y) & (x >= wx[0] + vertex_error["x"]) & (x <= wx[1] + vertex_error["x"]) & (y >= wy[0]) & (y <= wy[1])),
              "y_low": (base & np.isfinite(x) & np.isfinite(y) & (x >= wx[0]) & (x <= wx[1]) & (y >= wy[0] - vertex_error["y"]) & (y <= wy[1] - vertex_error["y"])),
              "y_high": (base & np.isfinite(x) & np.isfinite(y) & (x >= wx[0]) & (x <= wx[1]) & (y >= wy[0] + vertex_error["y"]) & (y <= wy[1] + vertex_error["y"]))
            }

            if cut_h["nominal"].sum() < 500:
                print(f"  {R} ohm {E:>4} GeV: solo {cut_h['nominal'].sum()} eventi dopo il taglio hodo, salto")
                continue

            out = dict(resistance=R, energy=E, energy_true=VERA.get(E, float(E)),
                       window="+".join(fallback) + "-ecal_prof" if fallback
                              else a.window,
                       n_centroid=int(cut_c["nominal"].sum()), n_hodo=int(cut_h["nominal"].sum()),
                       x_lo=wx[0], x_hi=wx[1], x_drop=wx[2], n_x=wx[3],
                       y_lo=wy[0], y_hi=wy[1], y_drop=wy[2], n_y=wy[3])

            for tag, cut_dict in (("cen", cut_c), ("hodo", cut_h)):

                variations_only = list(cut_dict.keys() - {"nominal"})
                variations_plus_nominal = variations_only + ["nominal"] #assures nominal is the last! Ruben

                reso_variations = []

                for var in variations_plus_nominal:

                  # code del DCB dal fit cumulativo di questa energia e di questa
                  # selezione: ad alta statistica sono determinate, run per run no
                  tails = None
                  if a.tails in ("fixed", "both") and cut_dict[var].sum() >= 500:
                      fp = fit_dcb(a_3x3[cut_dict[var]], E, R)
                      if fp is not None:
                          tails = fp["tails"]

                  def per_run(fix):
                      """{run: (sigma/mu, suo errore, sigma grezza, picco, suo errore,
                      eventi)} -- un dizionario, cosi' le due varianti del modello di fit
                      si allineano per run anche se una di loro perde qualche fit"""
                      o = {}
                      for r in sorted(int(u) for u in np.unique(run[cut_dict[var]])):
                          q = cut_dict[var] & (run == r)
                          if q.sum() < 300:
                              continue
                          fc = fit_dcb(a_3x3[q], E, R, fix=fix)
                          if fc is None:
                              continue

                          v, e = rel(fc)
                          o[r] = (v, e, rel(fc)[0], fc["peak"], fc["err_peak"],
                                  int(fc["nev"]))
                      return o

                  free = per_run(tails if a.tails == "fixed" else None)
                  fixed = per_run(tails) if (a.tails == "both" and tails) else {}
                  rns = sorted(free)
                  vals = [free[r][0] for r in rns]; errs = [free[r][1] for r in rns]
                  raws = [free[r][2] for r in rns]
                  peaks = [free[r][3] for r in rns]; epeaks = [free[r][4] for r in rns]
                  nevs = [free[r][5] for r in rns]
                  alt = [fixed[r][0] for r in rns if r in fixed]

                  # Dove nessun run singolo arriva a 300 eventi il punto non si perde:
                  # si fa UN fit cumulativo su tutti i run insieme. Succede a 340 ohm
                  # 250 GeV, dove il taglio sull'odoscopio lascia 846 eventi su 8 run.
                  pooled = False
                  fx = tails if a.tails == "fixed" else None
                  if not vals and cut_dict[var].sum() >= 500:
                      fc = fit_dcb(a_3x3[cut_dict[var]], E, R, fix=fx)
                      if fc is not None:
                          if fc is not None:
                              v, e = rel(fc)
                              vals, errs, raws, rns = [v], [e], [rel(fc)[0]], [0]
                              peaks, epeaks = [fc["peak"]], [fc["err_peak"]]
                              nevs = [int(fc["nev"])]
                              pooled = True
                  out[tag + "_pooled"] = int(pooled)
                  # pesi 1/sigma^2, come chiesto: l'errore statistico e' la sigma della
                  # sigma dalla varianza pesata con quei pesi
                  wts = [1.0 / (e * e) if e > 0 else 0.0 for e in errs]
                  mu, er = wmean(vals, errs, wts)
                  if var != "nominal": reso_variations.append(mu)

                reso_vtx_err_syst = np.sqrt( np.sum( (np.asarray(reso_variations) - mu)**2 )/4 ) #subtract nominal (- mu)  [Ruben]
                out[f"{tag}_vtx_syst"] = reso_vtx_err_syst
                print(R, E, "reso_vtx_err_syst: ", reso_vtx_err_syst, "nominal: ", mu)

                if not vals:
                    for k, v in ((tag, np.nan), (tag + "_stat", np.nan),
                                 (tag + "_drift", 0.),
                                 (tag + "_chi2", np.nan), (tag + "_tails", 0.),
                                 (tag + "_err", np.nan)):
                        out[k] = v
                    out["nrun_" + tag] = 0
                    out[tag + "_runs"] = []
                    continue


                sc = wscatter(vals, wts)

                est = max(er, sc) if np.isfinite(sc) else er
                out[tag] = mu
                out[tag + "_stat"] = est
                # sistematica sul modello di fit: quanto si sposta il punto se le code
                # del DCB si tengono fisse invece che libere. Non e' una larghezza da
                # togliere, e' un'ambiguita' sulla misura, quindi va nella barra
                ts = 0.
                if alt:
                    k = [r for r in rns if r in fixed]
                    ae = [free[r][1] for r in k]
                    aw = [1.0 / (e * e) if e > 0 else 0.0 for e in ae]
                    ma = wmean(alt, ae, aw)[0]
                    if np.isfinite(ma):
                        ts = abs(mu - ma)

                # DRIFT SUL PICCO, non sulla sigma. La sigma di ogni run e' misurata
                # attorno al picco di quel run, quindi il drift fra run non la sporca
                # direttamente: quello che si misura sui picchi e' l'instabilita' della
                # risposta, ed e' la stessa instabilita' che, agendo dentro un run,
                # allarga la sigma. Si stima dai picchi per run -- l'errore in piu' che
                # porta il loro fit a una costante a chi2/ndf = 1 -- e si esprime in
                # punti percentuali del picco medio, cosi' e' sottraibile in quadratura
                # da sigma/mu. Con un solo run non c'e' dispersione ed e' zero.
                dr, c0 = 0., np.nan
                if len(peaks) > 1:
                    qq = syst_for_unit_chi2(np.array(vals), np.array(errs))
                    c0 = float(qq[1])
                    dr = qq[0]
                out[tag + "_drift"] = dr
                out[tag + "_chi2"] = c0
                # la barra porta lo statistico e la sistematica sul modello di fit;
                # il drift NON sta nella barra, si sottrae
                out[tag + "_err"] = (dr**2 + est**2 + reso_vtx_err_syst**2)**0.5
                out["nrun_" + tag] = 0 if pooled else len(vals)
                out[tag + "_runs"] = [
                    dict(run=rn, nev=nv, sigma=v, err=e, peak=pk, err_peak=ep,
                         sigma_fixed=(fixed[rn][0] if rn in fixed else np.nan))
                    for rn, nv, v, e, pk, ep in zip(rns, nevs, vals, errs, peaks, epeaks)]
            b = bes_cons.get(E, 0.); syn = SYNC_C * out["energy_true"] ** 2.5
            print("bes_cons: ", b)
            b_nom = bes_nom.get(E, 0.)
            out["bes_cons"] = b; out["sync"] = syn; out["bes_nom"] = b_nom
            for tag in ("cen", "hodo"):
                # sottrazioni in quadratura: BES, sincrotrone, POS_eff (solo con il
                # centroide) e il drift dove esiste, cioe' dove il punto ha >1 run
                v = out[tag]; dr = out.get(tag + "_drift", 0.)
                out[tag + "_corr"] = (
                    math.sqrt(max(v * v - b * b - syn * syn, 0.))
                    if np.isfinite(v) else np.nan)
                out[tag + "_corr_larger_bes"] = (
                    math.sqrt(max(v * v - b_nom * b_nom - syn * syn, 0.))
                    if np.isfinite(v) else np.nan)
                out[tag + "_corr_larger_syn"] = (
                    math.sqrt(max(v * v - b * b - syn * syn * 1.3*1.3, 0.))
                    if np.isfinite(v) else np.nan)
                out[tag + "_bes_syst"] = abs(out[tag + "_corr_larger_bes"] - out[tag + "_corr"])
                out[tag + "_syn_syst"] = abs(out[tag + "_corr_larger_syn"] - out[tag + "_corr"])

                out[tag + "_err"] = (out[tag + "_err"]**2 + ( out[tag + "_bes_syst"] )**2)**0.5
                out[tag + "_err"] = (out[tag + "_err"]**2 + ( out[tag + "_syn_syst"] )**2)**0.5
            rows.append(out)
            print(f"  {R} ohm {E:>4} GeV: x [{wx[0]:+6.2f},{wx[1]:+6.2f}] mm "
                  f"(calo {wx[2]:.2f}%), y [{wy[0]:+6.2f},{wy[1]:+6.2f}] mm "
                  f"(calo {wy[2]:.2f}%) | centroide {out['cen']:.4f}% su "
                  f"{out['n_centroid']} ev, hodo {out['hodo']:.4f}% su {out['n_hodo']} ev "
                  f"({100*out['n_hodo']/max(out['n_centroid'],1):.0f}%)", flush=True)

    if not rows:
        return
    cols = ("resistance,energy,energy_true,window,n_centroid,n_hodo,nrun_cen,nrun_hodo,"
            "x_lo,x_hi,x_drop,n_x,y_lo,y_hi,y_drop,n_y,bes_cons,bes_nom,sync,"
            "cen,cen_stat,cen_drift,cen_chi2,cen_err,cen_corr,"
            "cen_pooled,hodo,hodo_stat,hodo_drift,hodo_chi2,"
            "hodo_err,hodo_corr,hodo_pooled,hodo_corr_larger_bes,hodo_bes_syst,hodo_vtx_syst")
    p = os.path.join(a.outdir, "resolution_hodo.csv")
    with open(p, "w") as fh:
        fh.write(cols + "\n")
        for r in rows:
            fh.write(",".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c])
                              for c in cols.split(",")) + "\n")
    print("->", p)

    # ---- tabella delle sistematiche, una riga per (resistenza, energia, catena) ----
    # Stessi numeri di resolution_hodo.csv ma in formato lungo e con il peso relativo
    # di ogni contributo, cosi' si legge a colpo d'occhio quale domina dove.
    scols = ("resistance", "energy", "energy_true", "chain", "window", "n_events",
             "n_run", "pooled", "sigma_raw", "stat", "syst_tails", "err_total",
             "bes", "sync", "drift", "chi2_peak", "sigma_corr",
             "bes_frac", "sync_frac", "drift_frac",
             "stat_frac", "tails_frac")
    p = os.path.join(a.outdir, "systematics.csv")
    with open(p, "w") as fh:
        fh.write(",".join(scols) + "\n")
        for r in rows:
            for tag in ("cen", "hodo"):
                v = r.get(tag, np.nan)
                if not np.isfinite(v) or v <= 0:
                    continue
                d = dict(resistance=r["resistance"], energy=r["energy"],
                         energy_true=r["energy_true"], chain=tag,
                         window=r.get("window", a.window),
                         n_events=r["n_centroid"] if tag == "cen" else r["n_hodo"],
                         n_run=r.get("nrun_" + tag, 0),
                         pooled=r.get(tag + "_pooled", 0),
                         sigma_raw=v, stat=r.get(tag + "_stat", np.nan),
                         syst_tails=r.get(tag + "_tails", 0.),
                         err_total=r.get(tag + "_err", np.nan),
                         bes=r["bes_cons"], sync=r["sync"],
                         drift=r.get(tag + "_drift", 0.),
                         chi2_peak=r.get(tag + "_chi2", np.nan),
                         sigma_corr=r.get(tag + "_corr", np.nan))
                # peso di ogni contributo in percentuale della sigma misurata
                for k, src in (("bes_frac", "bes"), ("sync_frac", "sync"),
                               ("drift_frac", "drift"),
                               ("stat_frac", "stat"), ("tails_frac", "syst_tails")):
                    d[k] = 100. * d[src] / v if np.isfinite(d[src]) else np.nan
                fh.write(",".join(f"{d[c]:.6g}" if isinstance(d[c], float) else str(d[c])
                                  for c in scols) + "\n")
    print("->", p)

    # ---- tutto per RUN ----------------------------------------------------
    # Il drift e' per energia per costruzione -- e' l'errore in piu' che rende
    # compatibili fra loro i picchi dei run di quell'energia -- ma cio' che lo genera
    # e' per run: lo scostamento del picco di ogni run dal picco medio. Qui ci sono
    # entrambi: le quantita' misurate run per run, e accanto i termini di energia a
    # cui contribuiscono, cosi' si vede quale run tira il drift.
    rcols = ("resistance", "energy", "energy_true", "chain", "run", "n_events",
             "sigma_pct", "err_sigma", "sigma_fixed_tails", "d_tails",
             "peak_ADC", "err_peak_ADC", "peak_dev_pct", "peak_pull",
             "sigma_dev_pct", "sigma_pull",
             "E_sigma_mean", "E_peak_mean_ADC", "E_drift_pct", "E_chi2_peak",
             "E_n_run", "E_bes", "E_sync", "E_stat", "E_syst_tails")
    p = os.path.join(a.outdir, "per_run.csv")
    with open(p, "w") as fh:
        fh.write(",".join(rcols) + "\n")
        for r in rows:
            for tag in ("cen", "hodo"):
                pr = r.get(tag + "_runs", [])
                if not pr:
                    continue
                pk = np.array([u["peak"] for u in pr])
                ep = np.array([u["err_peak"] for u in pr])
                w = 1. / np.maximum(ep, 1e-12) ** 2
                pm = float((pk * w).sum() / w.sum())
                sm = r.get(tag, np.nan)
                for u in pr:
                    d = dict(resistance=r["resistance"], energy=r["energy"],
                             energy_true=r["energy_true"], chain=tag, run=u["run"],
                             n_events=u["nev"], sigma_pct=u["sigma"],
                             err_sigma=u["err"], sigma_fixed_tails=u["sigma_fixed"],
                             d_tails=abs(u["sigma"] - u["sigma_fixed"]),
                             peak_ADC=u["peak"], err_peak_ADC=u["err_peak"],
                             peak_dev_pct=100. * (u["peak"] / pm - 1.),
                             peak_pull=(u["peak"] - pm) / u["err_peak"]
                                       if u["err_peak"] > 0 else np.nan,
                             sigma_dev_pct=100. * (u["sigma"] / sm - 1.)
                                           if sm > 0 else np.nan,
                             sigma_pull=(u["sigma"] - sm) / u["err"]
                                        if u["err"] > 0 else np.nan,
                             E_sigma_mean=sm, E_peak_mean_ADC=pm,
                             E_drift_pct=r.get(tag + "_drift", 0.),
                             E_chi2_peak=r.get(tag + "_chi2", np.nan),
                             E_n_run=r.get("nrun_" + tag, 0),
                             E_bes=r["bes_cons"], E_sync=r["sync"],
                             E_stat=r.get(tag + "_stat", np.nan),
                             E_syst_tails=r.get(tag + "_tails", 0.))
                    fh.write(",".join(f"{d[c]:.6g}" if isinstance(d[c], float)
                                      else str(d[c]) for c in rcols) + "\n")
    print("->", p)

    Rs = [R for R in (340, 400, 500) if any(r["resistance"] == R for r in rows)]

    # ---- controllo del drift, energia per energia -------------------------
    # Un pannello per energia: il PICCO dei singoli run, normalizzato al picco medio,
    # con il suo errore di fit, la costante e la banda di drift. Il drift e' l'errore
    # in piu' che porta il chi2/ndf di quel fit a 1: dove il chi2/ndf e' gia' <= 1 il
    # drift e' zero perche' non serve niente, e il pannello lo mostra.
    for tag, lab in (("cen", "centroid"), ("hodo", "hodoscope")):
        for R in Rs:
            q = [r for r in rows if r["resistance"] == R and len(r.get(tag + "_runs", [])) > 1]
            if not q:
                continue
            nc = 4; nr = int(np.ceil(len(q) / nc))
            fig, axs = plt.subplots(nr, nc, figsize=(4.4 * nc, 3.2 * nr), squeeze=False)
            for k, r in enumerate(q):
                ax = axs[k // nc][k % nc]
                pr = r[tag + "_runs"]
                xs = np.arange(len(pr))
                pk = np.array([u["peak"] for u in pr])
                ep = np.array([u["err_peak"] for u in pr])
                w = 1. / ep ** 2
                m0 = float((pk * w).sum() / w.sum())
                ys = 100 * (pk / m0 - 1.); es = 100 * ep / m0
                dr = r.get(tag + "_drift", 0.)
                c0 = r.get(tag + "_chi2", np.nan)
                ax.errorbar(xs, ys, yerr=es, fmt="o", ms=5, capsize=3, color="C0")
                ax.axhline(0., color="C3", lw=1.6,
                           label=f"weighted mean {m0:.1f} ADC")
                if dr > 0:
                    ax.axhspan(-dr, dr, color="C3", alpha=.15,
                               label=f"$\\pm$ drift {dr:.4f} %")
                ax.set_xticks(xs)
                ax.set_xticklabels([str(u["run"]) for u in pr], rotation=90, fontsize=6)
                ndf = len(pr) - 1
                if dr > 0:
                    txt = (f"$\\chi^2$/ndf {c0:.2f} ({ndf} ndf) $>$ 1  "
                           f"$\\rightarrow$ drift {dr:.4f} %")
                    col = "C3"
                else:
                    txt = (f"$\\chi^2$/ndf {c0:.2f} ({ndf} ndf) $\\leq$ 1  "
                           f"$\\rightarrow$ drift 0")
                    col = "C2"
                ax.set_title(f"{r['energy']} GeV\n{txt}", fontsize=8.5, color=col,
                             linespacing=1.5)
                ax.set_ylabel("peak / $\\langle$peak$\\rangle$ $-$ 1  [%]", fontsize=8)
                ax.tick_params(axis="y", labelsize=7); ax.grid(alpha=.3, axis="y")
                ax.legend(fontsize=6.5, loc="best")
            for k in range(len(q), nr * nc):
                axs[k // nc][k % nc].set_axis_off()
            fig.suptitle(f"{R} $\\Omega$ — cut on the {lab} — per-run PEAK against a "
                         f"constant.  green: already compatible, drift = 0.  "
                         f"red: extra error needed", fontsize=11)
            fig.tight_layout(rect=(0, 0, 1, 0.985))
            p = os.path.join(a.outdir, f"drift_check_{tag}_{R}ohm.png")
            fig.savefig(p, dpi=130); plt.close(fig)
            print("->", p)

    # due versioni dello stesso grafico: una con tutti i punti dentro al fit, una in
    # cui --nofit-energies restano disegnati ma fuori dal fit. Cosi' si vede subito
    # quanto pesano quei punti sui parametri
    for tag, lab in (("cen", "centroid"), ("hodo", "hodoscope")):
      for allfit in (False, True):
          fig, axs = plt.subplots(2, len(Rs), figsize=(6.4 * len(Rs), 10), sharex=True,
                                  gridspec_kw=dict(height_ratios=[2, 1]), squeeze=False)
          S_fixed, C_fixed = None, None
          for j, R in enumerate(Rs):
              q = [r for r in rows if r["resistance"] == R]
              x = np.array([r["energy_true"] for r in q])
              raw = np.array([r[tag] for r in q])
              y = np.array([r[tag + "_corr"] for r in q])
              e = np.array([r[tag + "_err"] for r in q])
              bes = np.array([r["bes_cons"] for r in q]); syn = np.array([r["sync"] for r in q])
              dri = np.array([r.get(tag + "_drift", 0.) for r in q])
              bes_syst = np.array([r[tag + "_bes_syst"] for r in q])
              syn_syst = np.array([r[tag + "_syn_syst"] for r in q])
              vtx_syst = np.array([r[tag + "_vtx_syst"] for r in q])
              ax, ax2 = axs[0][j], axs[1][j]
              ax.plot(x, raw, "o-", ms=6, color="0.35", label="$\\sigma/\\mu$")
              ax.errorbar(x, y, yerr=e, fmt="^-", ms=7.5, color="C3", capsize=3,
                          alpha=.9, label="$-$ BES $-$ sync" )
              infit = (np.ones(len(q), bool) if allfit else
                       ~np.isin([r["energy"] for r in q], a.nofit_energies))
              g = np.isfinite(y) & (y > 0) & (e > 0) & infit
              shown = np.isfinite(y) & (y > 0) & (e > 0) & ~infit
              if shown.any():
                  ax.errorbar(x[shown], y[shown], yerr=e[shown], fmt="^", ms=7.5,
                              mfc="none", color="C3", capsize=3, label="not in the fit")
              pl = np.array([bool(r.get(tag + "_pooled", 0)) for r in q])
              pl &= np.isfinite(y) & (y > 0)
              if pl.any():
                  ax.plot(x[pl], y[pl], "s", ms=14, mfc="none", mew=1.4, color="C4",
                          label="one pooled fit (no run has enough events)")
              fb = np.array([str(r.get("window", "")).endswith("-ecal_prof") for r in q])
              fb &= np.isfinite(y) & (y > 0)
              if fb.any() and tag == "hodo":
                  ax.plot(x[fb], y[fb], "o", ms=13, mfc="none", mew=1.4, color="0.25",
                          label="no parabola: used ECAL/hodo profile")
              if g.sum() >= 4:
                  mi = Minuit(LeastSquares(x[g], y[g], e[g], reso), N=0.3, S=3., C=0.3)
                  for k in "NSC":
                      mi.limits[k] = (0, None)
                  if R != 340 and (C_fixed is not None) and (S_fixed is not None):
                      mi.values["S"] = S_fixed; mi.fixed["S"] = True
                      #mi.values["C"] = C_fixed; mi.fixed["C"] = True
                  mi.migrad(); mi.hesse()
                  if R == 340:
                    S_fixed = mi.values['S']
                    C_fixed = mi.values['C']
                  xs = np.linspace(x[g].min() * .9, x[g].max() * 1.05, 300)
                  ax.plot(xs, reso(xs, *mi.values), "--", lw=2.2, color="darkviolet",
                          label="fit  $N/E \\oplus S/\\sqrt{E} \\oplus C$")
                  nd = int(g.sum()) - (2 if R == 500 else 3)
                  ax.text(.77, .75,
                          f"$N$ {1000*mi.values['N']:6.0f} $\\pm$ {1000*mi.errors['N']:.0f} MeV\n"
                          f"$S$ {mi.values['S']:6.3f}" +
                          (" % (FIXED)" if R != 340 else f" $\\pm$ {mi.errors['S']:.3f} %") +
                          f"\n$C$ {mi.values['C']:6.3f}" +
                          (" % (FIXED)" if R != 340 else f" $\\pm$ {mi.errors['C']:.3f} %") +
                          f"\n$\\chi^2$/ndf {mi.fval:6.1f} / {nd}",
                          transform=ax.transAxes, ha="right", va="top", fontsize=9.5,
                          family="monospace", bbox=dict(fc="white", ec="darkviolet", pad=6))
              ax.set_title(f"{R} $\\Omega$", fontsize=12, fontweight="bold")
              ax.set_ylabel("$\\sigma/E$  [%]"); ax.grid(alpha=.3); ax.legend(fontsize=8)
              if np.isfinite(np.nanmax(raw)):
                  ax.set_ylim(0, float(np.nanmax(raw) * 1.08))
              ax2.plot(x, raw, "o-", ms=5, color="0.35", label="$\\sigma/\\mu$")
              ax2.plot(x, bes, "D-.", ms=5, color="C1", label="BES")
              ax2.plot(x, syn, "^-", ms=5, color="C4", label="synchrotron")
              nr = np.array([r.get("nrun_" + tag, 0) for r in q])
              ax2.plot(x, np.where(dri > 0, dri, np.nan), "s--", ms=5, color="C0",
                       label="drift")
              ax2.plot(x, np.where(vtx_syst > 0, vtx_syst, np.nan), "s--", ms=5, color="C5",
                       label="vtx syst")
              ax2.plot(x, np.where(bes_syst > 0, bes_syst, np.nan), "s--", ms=5, color="C6",
                       label="BES syst")
              ax2.plot(x, np.where(syn_syst > 0, syn_syst, np.nan), "s--", ms=5, color="C7",
                       label="Synchr. syst")
              # dove il drift non c'e' si dice perche': un solo run (niente dispersione
              # da misurare) oppure sigma per run gia' compatibili fra loro
              lo2 = float(np.nanmin(np.concatenate([bes, syn, [1e-3]])))
              m1 = (dri <= 0) & (nr <= 1)
              mc = (dri <= 0) & (nr > 1)
              if m1.any():
                  ax2.plot(x[m1], np.full(m1.sum(), lo2), "x", ms=6, color="C0",
                           label="drift n/a (1 run)")
              if mc.any():
                  ax2.plot(x[mc], np.full(mc.sum(), lo2), "s", ms=5, mfc="none",
                           color="C0", label="drift = 0 ($\\chi^2$/ndf $\\leq$ 1)")
              t = np.concatenate([raw, bes, syn, dri[dri > 0], [1e-3]])
              t = t[np.isfinite(t) & (t > 0)]
              ax2.set_yscale("log")
              if len(t):
                  ax2.set_ylim(float(t.min() * .4), float(t.max() * 2.5))
              ax2.set_xlabel("True beam energy [GeV]")
              ax2.set_ylabel("size of each term  [%]")
              ax2.grid(alpha=.3, which="both"); ax2.legend(fontsize=8)
          fig.suptitle(f"cut on the {lab}   $\\quad$   $A_{{tot}} > {A_TOT_MIN:.0f}$ ADC"
                       + ("   $\\quad$   all points in the fit" if allfit else
                          "   $\\quad$   "),  fontsize=12)

          fig.tight_layout()
          p = os.path.join(a.outdir,
                           f"resolution_terms_{tag}{'_allpoints' if allfit else ''}.png")
          fig.savefig(p, dpi=150); plt.close(fig)
          print("->", p)


if __name__ == "__main__":
    main()

