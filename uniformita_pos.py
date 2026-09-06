"""
Response non-uniformity in (pos_eta, pos_phi): four treatments compared.

The response A_tot(eta, phi) inside the selection window is not flat: it is a
surface with negative curvature in both coordinates, more pronounced in eta than in
phi, and the curvature decreases with energy. It is measured by profili_pernorm.py,
which writes c/a in % per crystal^2 for every (resistance, energy).
The beam illuminates it non-uniformly, and differently from run to run, so part of
the measured width comes from the impact position and not from the calorimeter.

Note: making the illumination FLAT does not remove the term. With uniform
illumination over |u| <= 0.2 and relative curvature c, the residual spread of the
response factor is c * std(u^2) = c * 0.0119. Flat reweighting standardises the
reference, it does not correct. To remove the term the event has to be corrected.

Treatments (all in the runmean chain: DCB fit per run, then the mean of
sigma_i / peak_i weighted by the number of events):

  raw    no treatment, reference
  corr   A -> A * <f>_run / f(u_i, v_i), event-by-event correction with the response
         surface. Removes the term and loses no statistics.
  pos    no correction to the fit; the term is subtracted in quadrature, per run,
         with POS_r = 100 * std(f_i)/mean(f_i) over the events of that run:
         sigma_corr = sqrt((sigma/mu)^2 - POS^2)
  flat   reweighting to flat illumination, w proportional to 1/occupancy on a
         --grid x --grid mesh, weighted DCB fit per run

Response surface f: full quadratic f = a0 + a1 u + a2 v + a3 u^2 + a4 v^2 + a5 uv
(u = pos_eta - 18, v = pos_phi - 6), by LINEAR least squares on the binned 2D map ->
unique solution, no vertex/curvature degeneracy (a shift of the vertex is
compensated by a change of curvature). Estimated PER ENERGY on events normalised to
the peak of their own run, so that the two response populations at 340 ohm do not
enter the map.

Usage:
  python3 plot/uniformita_pos.py --base . --outdir plot/uniformita \
      --resistances 340 --energies 20 30 40 --exclude-runs 20592 ...
  python3 plot/uniformita_pos.py --only-collect --besdir plot/bes
"""

import argparse, os, glob, re, json, math
import runsets
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares

SCALE = {340: 3500/150., 400: 1080/40., 500: 3340/100.}
NBINS, XLO, XHI = 8000, 0., 8000.
ETA0, PHI0, SEL = 18., 6., 0.2
A_TOT_MIN = 100.
# the map and the POS term are built on the same events that enter the sigma fit,
# i.e. inside the fit.sh window of the run: using +-10 sigma inflates POS, because the
# tails sit at the edges in (u, v).
NMIN_BIN = 1              # eventi minimi per bin della mappa
SYNC_C = 1.92e-7

CUT = ("$|\\mathrm{pos\\_eta}-18| \\leq 0.2$,  $|\\mathrm{pos\\_phi}-6| \\leq 0.2$")

VERA = {20:20.00, 30:30.00, 40:39.99, 50:49.98, 60:59.97, 80:79.90, 100:99.75,
        120:119.48, 150:148.73, 175:172.67, 200:196.08, 225:218.82, 250:240.76,
        275:261.77, 300:281.74}

FILES = {340: ("reco_340ohm", "*_merged.root"),
         400: ("reco_400ohm", "*_400_merged.root"),
         500: ("reco_500ohm", "*_500_merged.root")}

PARS = ("alpha_l", "alpha_h", "n_l", "n_h", "mean", "sigma", "N")

METHODS = ("raw", "corr", "pos", "flat")
MCOL = {"raw": "k", "corr": "C3", "pos": "C0", "flat": "C2"}
MMK = {"raw": "o", "corr": "^", "pos": "s", "flat": "v"}


# ------------------------------------------------------------------- DCB
def dcb_func(x, alpha_l, alpha_h, n_l, n_h, mean, sigma, N):
    t = (x - mean) / sigma
    out = np.empty_like(t)
    core = (t >= -alpha_l) & (t <= alpha_h)
    low, high = t < -alpha_l, t > alpha_h
    out[core] = np.exp(-0.5 * t[core] ** 2)
    if np.any(low):
        f2 = (n_l/alpha_l) - alpha_l - t[low]
        out[low] = np.exp(-0.5*alpha_l**2) * np.power(np.maximum(alpha_l/n_l*f2, 1e-12), -n_l)
    if np.any(high):
        f2 = (n_h/alpha_h) - alpha_h + t[high]
        out[high] = np.exp(-0.5*alpha_h**2) * np.power(np.maximum(alpha_h/n_h*f2, 1e-12), -n_h)
    return N * out


def hist_stats(c, x, lo, hi):
    m = (x >= lo) & (x <= hi)
    cc, xx = c[m], x[m]
    tot = cc.sum()
    if tot <= 0:
        return np.nan, np.nan, 0.
    mu = (cc*xx).sum()/tot
    return mu, np.sqrt(max((cc*(xx-mu)**2).sum()/tot, 0.)), tot


def fit_window(v, energy, resistance):
    c, e = np.histogram(v, bins=NBINS, range=(XLO, XHI))
    c = c.astype(float); x = 0.5*(e[:-1]+e[1:])
    sc = SCALE[resistance]
    lo, hi = sc*energy*0.95, sc*energy*1.05
    for _ in range(2):
        mu, rms, tot = hist_stats(c, x, lo, hi)
        if not np.isfinite(mu) or rms <= 0:
            return None
        lo, hi = mu-3*rms, mu+3*rms
    return lo, hi


def mode_window(v, energy, resistance):
    nom = SCALE[resistance]*energy
    x = v[(v > 0.5*nom) & (v < 1.3*nom)]
    if len(x) < 100:
        return None
    c, e = np.histogram(x, bins=150)
    mode = 0.5*(e[c.argmax()] + e[c.argmax()+1])
    core = x[np.abs(x-mode) < 0.08*nom]
    if len(core) < 50 or core.std() <= 0:
        return None
    return mode-3*core.std(), mode+3*core.std()


def _win_ok(w, v):
    if w is None:
        return False
    lo, hi = w
    return hi > lo and ((v >= lo) & (v <= hi)).sum() >= 200


def fd_binwidth(v, n_eff=None):
    if len(v) < 20:
        return 1.
    q25, q75 = np.percentile(v, [25, 75])
    n = n_eff if n_eff else len(v)
    return float(max(1., round(2*(q75-q25)/max(n, 1)**(1./3.))))


TAILS = ("alpha_l", "alpha_h", "n_l", "n_h")


def fit_dcb(v, energy, resistance, w=None, fix=None):
    """DCB fit, optionally weighted. The fit.sh window, with Freedman-Diaconis binning
    ALWAYS computed on the unweighted events, so that the weighted and unweighted
    versions have the same ndf and can be compared.

    fix: dict of tail parameters (alpha_l, alpha_h, n_l, n_h) to hold constant. Left
    free they are badly determined -- n_l and n_h rail against their limit of 10 in
    most runs -- and sigma, being correlated with them, inherits an error about twice
    what it would otherwise have. Holding them at the values of the pooled fit of the
    same energy shrinks the error on sigma by 1.4 to 1.9, and a bootstrap confirms the
    smaller error is the true one. It is not free: sigma itself moves by up to 1 %, so
    it is a change of fit model and is carried as a systematic, not adopted silently."""
    win = fit_window(v, energy, resistance)

    if not _win_ok(win, v):
        win = mode_window(v, energy, resistance)
    if not _win_ok(win, v):
        return None
    lo, hi = win
    inwin = v[(v >= lo) & (v <= hi)]
    if len(inwin) < 200:
        return None
    bw = fd_binwidth(inwin)
    nb = max(int(round((hi-lo)/bw)), 12)
    edges = np.linspace(lo, hi, nb+1)
    x = 0.5*(edges[:-1]+edges[1:])
    if w is None:
        y, _ = np.histogram(v, bins=edges)
        y = y.astype(float)
        ey = np.sqrt(np.maximum(y, 1.))
    else:
        y, _ = np.histogram(v, bins=edges, weights=w)
        y2, _ = np.histogram(v, bins=edges, weights=w**2)
        ey = np.sqrt(np.maximum(y2, 0.))
    sel = (y > 0) & np.isfinite(ey) & (ey > 0)
    if sel.sum() < 12:
        return None
    x, y, ey = x[sel], y[sel], ey[sel]
    seed = dict(alpha_l=2., alpha_h=2., n_l=2., n_h=2.,
                mean=(y*x).sum()/y.sum(), sigma=0.5*(hi-lo)/3., N=float(y.max()))
    if fix:
        seed.update({k: fix[k] for k in TAILS if k in fix})
    best = None
    for _ in range(3):
        m = Minuit(LeastSquares(x, y, ey, dcb_func), **seed)
        m.limits["alpha_l"] = (0.1, 10); m.limits["alpha_h"] = (0.1, 10)
        m.limits["n_l"] = (1, 10); m.limits["n_h"] = (1, 10)
        m.limits["mean"] = (lo, hi); m.limits["sigma"] = (0, hi-lo); m.limits["N"] = (0, None)
        if fix:
            for k in TAILS:
                if k in fix:
                    m.fixed[k] = True
        m.migrad(); m.hesse(); best = m
        seed = {p: m.values[p] for p in seed}
    if (best.covariance is None or best.errors["sigma"] <= 0 or best.errors["mean"] <= 0):
        m = Minuit(LeastSquares(x, y, ey, dcb_func), **{p: best.values[p] for p in PARS})
        for p in ("alpha_l", "alpha_h", "n_l", "n_h"):
            m.limits[p] = best.limits[p]; m.fixed[p] = True
        m.limits["mean"] = (lo, hi); m.limits["sigma"] = (0, hi-lo); m.limits["N"] = (0, None)
        m.migrad(); m.hesse()
        if m.errors["sigma"] > 0 and m.errors["mean"] > 0:
            best = m
    return dict(peak=float(best.values["mean"]), err_peak=float(best.errors["mean"]),
                sigma=float(best.values["sigma"]), err_sigma=float(best.errors["sigma"]),
                chi2=float(best.fval), nev=int(len(inwin)),
                ndf=max(len(x) - (3 if fix else 7), 1),
                tails={k: float(best.values[k]) for k in TAILS},
                lo=float(lo), hi=float(hi))


def rel(r):
    v = 100*r["sigma"]/r["peak"]
    e = v*math.sqrt((r["err_sigma"]/r["sigma"])**2 + (r["err_peak"]/r["peak"])**2)
    return v, e


def wscatter(vals, wts):
    """Error on the weighted mean from the WEIGHTED VARIANCE of the values, not from the
    errors of the individual fits: it already contains both the fit noise and the
    run-to-run spread. n_eff = (sum w)^2 / sum w^2. With a single run it is
    undefined -> nan."""
    v = np.asarray(vals, float); w = np.asarray(wts, float)
    g = np.isfinite(v) & (w > 0)
    v, w = v[g], w[g]
    if len(v) < 2:
        return np.nan
    V1, V2 = w.sum(), (w ** 2).sum()
    neff = V1 ** 2 / V2
    if neff <= 1:
        return np.nan
    mu = (w * v).sum() / V1
    return float(np.sqrt((w * (v - mu) ** 2).sum() / V1 / (neff - 1)))


def wmean(vals, errs, wts):
    vals = np.asarray(vals, float); errs = np.asarray(errs, float); wts = np.asarray(wts, float)
    g = np.isfinite(vals) & np.isfinite(errs) & (wts > 0)
    if g.sum() == 0:
        return np.nan, np.nan
    v, e, w = vals[g], errs[g], wts[g]
    return float((v*w).sum()/w.sum()), float(np.sqrt((w**2*e**2).sum())/w.sum())


# ------------------------------------------------------- superficie di risposta
def design(u, v):
    return np.column_stack([np.ones_like(u), u, v, u*u, v*v, u*v])


def fit_surface(u, v, a, ngrid, nmin=NMIN_BIN):
    """Linear least squares of f = a0 + a1 u + a2 v + a3 u^2 + a4 v^2 + a5 uv on the
    binned 2D map of the means. Returns (coef, cov, diagnostics).

    Error per bin: sigma_pooled / sqrt(N_bin), with sigma_pooled = RMS of a inside
    the box. NOT the RMS of the individual bin: with fine meshes (150x150 over +-0.2
    is 22500 bins for ~30k events, i.e. ~1 event per bin) the per-bin RMS is not
    defined. The spread of a is the same everywhere in the box to within ~0.2%, so
    the pooled estimate is the right one, and at fine binning this fit coincides
    with the unbinned linear regression on the individual events."""
    rng = [[-SEL, SEL], [-SEL, SEL]]
    H, ue, ve = np.histogram2d(u, v, bins=ngrid, range=rng)
    S, _, _ = np.histogram2d(u, v, bins=ngrid, range=rng, weights=a)
    pooled = float(np.std(a))
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = S/np.maximum(H, 1)
        err = pooled/np.sqrt(np.maximum(H, 1))
    uc = 0.5*(ue[:-1]+ue[1:]); vc = 0.5*(ve[:-1]+ve[1:])
    UU, VV = np.meshgrid(uc, vc, indexing="ij")
    ok = (H >= nmin) & (err > 0)
    if ok.sum() < 12:
        return None
    X = design(UU[ok], VV[ok])
    y = mean[ok]; s = err[ok]
    A = X/s[:, None]; b = y/s
    ATA = A.T @ A
    coef = np.linalg.solve(ATA, A.T @ b)
    cov = np.linalg.inv(ATA)
    resid = (y - X @ coef)/s
    diag = dict(nbin=int(ok.sum()), chi2=float((resid**2).sum()),
                ndf=int(ok.sum()-6), nmin=float(H[ok].min()), nmax=float(H[ok].max()),
                pooled=pooled, nev_map=int(len(a)))
    return coef, cov, diag


def surf_eval(coef, u, v):
    return design(u, v) @ coef


# -------------------------------------------------------------------- pesi
def flat_weights(u, v, ngrid):
    rng = [[-SEL, SEL], [-SEL, SEL]]
    H, ue, ve = np.histogram2d(u, v, bins=ngrid, range=rng)
    iu = np.clip(np.digitize(u, ue)-1, 0, ngrid-1)
    iv = np.clip(np.digitize(v, ve)-1, 0, ngrid-1)
    occ = H[iu, iv]
    w = np.where(occ > 0, 1./np.maximum(occ, 1), 0.)
    if w.sum() <= 0:
        return None
    w = w*len(u)/w.sum()
    return w


# ---------------------------------------------------------------- analisi
def analyse(path, energy, resistance, drop, ngrid, gridflat, nmin=NMIN_BIN, only=()):
    t = uproot.open(path)["h4_reco"]
    arr = t.arrays(["run", "A_tot", "pos_eta", "pos_phi"], library="np")
    k = ((np.abs(arr["pos_eta"]-ETA0) <= SEL) & (np.abs(arr["pos_phi"]-PHI0) <= SEL)
         & (arr["A_tot"] > A_TOT_MIN))
    if drop:
        k &= ~np.isin(arr["run"], drop)
    if len(only):
        k &= np.isin(arr["run"], only)
    run = arr["run"][k]
    at = arr["A_tot"][k].astype(float)
    u = arr["pos_eta"][k].astype(float) - ETA0
    v = arr["pos_phi"][k].astype(float) - PHI0
    if len(at) < 500:
        return None
    runs = sorted(int(r) for r in np.unique(run))

    # --- passo 1: fit DCB per run, riferimento raw
    base = {}
    for r in runs:
        m = run == r
        f = fit_dcb(at[m], energy, resistance)
        if f is not None:
            base[r] = f

    if not base:
        return None

    # --- passo 2: mappa di risposta per energia, su eventi normalizzati al run
    mmap = np.zeros(len(at), bool)
    anorm = np.zeros(len(at))
    for r, f in base.items():
        m = run == r
        anorm[m] = at[m]/f["peak"]
        mmap[m] = (at[m] >= f["lo"]) & (at[m] <= f["hi"])
    surf = fit_surface(u[mmap], v[mmap], anorm[mmap], ngrid, nmin)
    if surf is None:
        return None
    coef, cov, sdiag = surf
    fev = surf_eval(coef, u, v)

    # relative curvature (after the normalisation: a0 ~ 1)
    a0 = float(coef[0])
    sdiag.update(coef=[float(c) for c in coef],
                 err_coef=[float(np.sqrt(cov[i, i])) for i in range(6)],
                 curv_eta_pct=100*float(coef[3])/a0,
                 curv_phi_pct=100*float(coef[4])/a0,
                 cross_pct=100*float(coef[5])/a0)

    # --- passo 3: i tre trattamenti, per run
    per = []
    for r in runs:
        m = run == r
        row = dict(run=r, n=int(m.sum()))
        b = base.get(r)
        if b is None:
            per.append(row); continue
        row["raw"] = rel(b); row["peak"] = b["peak"]
        row["chi2ndf_raw"] = b["chi2"]/b["ndf"]

        win = mmap[m]
        row["n_win"] = int(win.sum())
        fr = fev[m][win]
        fm = float(fr.mean())
        if fm > 0 and win.sum() > 100:
            row["pos_pct"] = 100*float(fr.std())/fm
            aw = at[m][win]
            row["rms_raw"] = 100*float(aw.std())/float(aw.mean())
            ac = aw*fm/np.maximum(fr, 1e-9)
            row["rms_corr"] = 100*float(ac.std())/float(ac.mean())
            row["rms_att"] = float(np.sqrt(max(row["rms_raw"]**2 - row["pos_pct"]**2, 0.)))
            fc = fit_dcb(at[m]*fm/np.maximum(fev[m], 1e-9), energy, resistance)
            if fc is not None:
                row["corr"] = rel(fc); row["chi2ndf_corr"] = fc["chi2"]/fc["ndf"]
        # flat
        wf = flat_weights(u[m], v[m], gridflat)
        if wf is not None:
            row["neff"] = float(wf.sum()**2/(wf**2).sum())
            ff = fit_dcb(at[m], energy, resistance, w=wf)
            if ff is not None:
                row["flat"] = rel(ff); row["chi2ndf_flat"] = ff["chi2"]/ff["ndf"]
        per.append(row)

    out = dict(energy=energy, resistance=resistance, nrun=len(runs), nev=int(len(at)),
               runs=runs, surface=sdiag, per_run=per, grid=ngrid, grid_flat=gridflat,
               nmin=nmin)

    ok = [p for p in per if "raw" in p]
    wts = [p["n"] for p in ok]
    out["raw"] = wmean([p["raw"][0] for p in ok], [p["raw"][1] for p in ok], wts)
    out["raw_scat"] = wscatter([p["raw"][0] for p in ok], wts)
    for meth in ("corr", "flat"):
        sub = [p for p in ok if meth in p]
        out[meth] = wmean([p[meth][0] for p in sub], [p[meth][1] for p in sub],
                          [p["n"] for p in sub]) if sub else (np.nan, np.nan)
        out[meth + "_scat"] = (wscatter([p[meth][0] for p in sub], [p["n"] for p in sub])
                               if sub else np.nan)
    # pos: sottrazione in quadratura per run, poi media
    sub = [p for p in ok if "pos_pct" in p]
    vals, errs, wp = [], [], []
    for p in sub:
        s, es = p["raw"]
        q = s*s - p["pos_pct"]**2
        if q <= 0:
            continue
        sc = math.sqrt(q)
        vals.append(sc); errs.append(es*s/sc); wp.append(p["n"])
    out["pos"] = wmean(vals, errs, wp) if vals else (np.nan, np.nan)
    out["pos_term"] = wmean([p["pos_pct"] for p in sub],
                            [0. for p in sub], [p["n"] for p in sub])[0] if sub else np.nan
    out["neff_frac"] = (float(np.sum([p.get("neff", 0.) for p in ok]) /
                              max(np.sum(wts), 1)) if ok else np.nan)
    for key in ("rms_raw", "rms_corr", "rms_att"):
        s2 = [p for p in ok if key in p]
        out[key] = (wmean([p[key] for p in s2], [0.]*len(s2), [p["n_win"] for p in s2])[0]
                    if s2 else np.nan)
    return out


# --------------------------------------------------------------- N,S,C
def load_bes(besdir, R):
    b = {}
    f = os.path.join(besdir, f"rereco_{R}_withBES.csv")
    if not os.path.exists(f):
        return b
    for i, l in enumerate(open(f)):
        if i == 0:
            continue
        p = l.strip().split(",")
        if len(p) >= 7:
            b[int(float(p[0]))] = float(p[6])
    return b


def nsc(x, N, S, C):
    return np.sqrt((100*N/x)**2 + (S/np.sqrt(x))**2 + C**2)


def fit_nsc(E, y, ey):
    m = Minuit(LeastSquares(E, y, ey, nsc), N=0.3, S=3., C=0.3)
    m.limits["N"] = (0, 2); m.limits["S"] = (0, 20); m.limits["C"] = (0, 5)
    m.migrad(); m.hesse()
    return m


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=".")
    ap.add_argument("--outdir", default="plot/uniformita")
    ap.add_argument("--besdir", default="plot/bes")
    ap.add_argument("--resistances", nargs="+", type=int, default=[340, 400, 500])
    ap.add_argument("--energies", nargs="*", type=int, default=None)
    ap.add_argument("--exclude-runs", nargs="*", type=int, default=[])
    ap.add_argument("--grid", type=int, default=150, help="griglia per la superficie f")
    ap.add_argument("--nmin", type=int, default=NMIN_BIN,
                    help="eventi minimi per bin della mappa")
    ap.add_argument("--grid-flat", type=int, default=10, help="griglia per i pesi piatti")
    runsets.add_argument(ap)
    ap.add_argument("--only-collect", action="store_true")
    a = ap.parse_args()
    drop, only = runsets.resolve(a.runset, a.exclude_runs)
    os.makedirs(a.outdir, exist_ok=True)
    cache = os.path.join(a.outdir, "_cache"); os.makedirs(cache, exist_ok=True)

    if not a.only_collect:
        for R in a.resistances:
            d, pat = FILES[R]
            files = sorted(glob.glob(os.path.join(a.base, d, pat)),
                           key=lambda p: int(re.match(r"(\d+)", os.path.basename(p)).group(1)))
            for f in files:
                E = int(re.match(r"(\d+)", os.path.basename(f)).group(1))
                if a.energies is not None and E not in a.energies:
                    continue
                print(f"  {R} ohm {E:4d} GeV", flush=True)
                res = analyse(f, E, R, drop, a.grid, a.grid_flat, a.nmin, only)
                if res is None:
                    print("      salto"); continue
                json.dump(res, open(os.path.join(cache, f"{R}_{E}.json"), "w"), default=float)
                s = res["surface"]
                print(f"      raw {res['raw'][0]:.4f}  corr {res['corr'][0]:.4f}  "
                      f"pos {res['pos'][0]:.4f} (POS {res['pos_term']:.4f})  "
                      f"flat {res['flat'][0]:.4f}   |   superficie "
                      f"c_eta {s['curv_eta_pct']:+.2f} c_phi {s['curv_phi_pct']:+.2f} "
                      f"c_x {s['cross_pct']:+.2f} %/cristallo^2, "
                      f"chi2/ndf {s['chi2']/max(s['ndf'],1):.2f}\n"
                      f"          RMS troncato: raw {res['rms_raw']:.4f} -> corr "
                      f"{res['rms_corr']:.4f} (quadratura {res['rms_att']:.4f})", flush=True)

    # ------------------------------------------------------------ raccolta
    rows = []
    for R in a.resistances:
        bes = load_bes(a.besdir, R)
        for cf in sorted(glob.glob(os.path.join(cache, f"{R}_*.json")),
                         key=lambda p: int(os.path.basename(p).split("_")[1].split(".")[0])):
            c = json.load(open(cf))
            E = int(c["energy"]); Et = VERA.get(E, float(E))
            syn = SYNC_C*Et**2.5
            b = bes.get(E, 0.)
            s = c["surface"]
            r = dict(resistance=R, energy=E, energy_true=Et, nrun=c["nrun"], nev=c["nev"],
                     bes=b, sync=syn, pos_term=c["pos_term"],
                     curv_eta=s["curv_eta_pct"], curv_phi=s["curv_phi_pct"],
                     cross=s["cross_pct"], surf_chi2ndf=s["chi2"]/max(s["ndf"], 1),
                     neff_frac=100*c.get("neff_frac", np.nan),
                     rms_raw=c.get("rms_raw", np.nan), rms_corr=c.get("rms_corr", np.nan),
                     rms_att=c.get("rms_att", np.nan),
                     raw_scat=c.get("raw_scat", np.nan),
                     corr_scat=c.get("corr_scat", np.nan))
            for m in METHODS:
                v, e = c[m]
                r[m] = v; r[m+"_err"] = e
                r[m+"_corr"] = (float(np.sqrt(max(v*v - b*b - syn*syn, 0.)))
                                if np.isfinite(v) else np.nan)
            trio = [r[m+"_corr"] for m in ("corr", "pos", "flat")]
            base = r["raw_corr"]
            dev = [abs(t-base) for t in trio if np.isfinite(t)]
            r["syst_pct"] = float(max(dev)) if dev else np.nan
            r["syst_rel"] = 100*r["syst_pct"]/base if (dev and base > 0) else np.nan
            rows.append(r)

    cols = (["resistance","energy","energy_true","nrun","nev","bes","sync","pos_term",
             "curv_eta","curv_phi","cross","surf_chi2ndf","neff_frac"]
            + [c for m in METHODS for c in (m, m+"_err", m+"_corr")]
            + ["syst_pct","syst_rel"])
    csv = os.path.join(a.outdir, "uniformita_pos.csv")
    with open(csv, "w") as fh:
        fh.write(",".join(cols)+"\n")
        for r in rows:
            fh.write(",".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c])
                              for c in cols)+"\n")
    print("->", csv)

    # --------------------------------------------------------------- plot
    fitres = {}
    for R in a.resistances:
        rr = [r for r in rows if r["resistance"] == R]
        if len(rr) < 2:
            continue
        fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True,
                                gridspec_kw=dict(height_ratios=[2.2, 1, 1]))
        box = []
        for m in METHODS:
            E = np.array([r["energy_true"] for r in rr])
            y = np.array([r[m+"_corr"] for r in rr])
            ey = np.array([r[m+"_err"] for r in rr])
            g = np.isfinite(y) & (y > 0) & np.isfinite(ey) & (ey > 0)
            axs[0].errorbar(E[g], y[g], yerr=ey[g], fmt=MMK[m]+"-", color=MCOL[m],
                            capsize=3, ms=6, lw=1.1, label=m)
            if g.sum() >= 4:
                mi = fit_nsc(E[g], y[g], ey[g])
                fitres[(R, m)] = (mi, int(g.sum()))
                box.append(f"{m:<5s} N={1000*mi.values['N']:5.0f} MeV  "
                           f"S={mi.values['S']:5.2f}%  C={mi.values['C']:6.3f}%  "
                           f"$\\chi^2$/ndf={mi.fval:6.1f}/{g.sum()-3}")
                xs = np.linspace(E[g].min(), E[g].max(), 300)
                axs[0].plot(xs, nsc(xs, *mi.values), "--", color=MCOL[m], lw=1, alpha=.6)
        axs[0].set_ylabel("$\\sigma/E$ [%]   ($\\ominus$ BES $\\ominus$ sync)")
        axs[0].set_title(f"{R} $\\Omega$   $\\quad$   $A_{{tot}} > {A_TOT_MIN:.0f}$ ADC,  {CUT}",
                         fontsize=12)
        axs[0].legend(fontsize=10, loc="lower left"); axs[0].grid(alpha=.3)
        if box:
            axs[0].text(.98, .95, "\n".join(box), transform=axs[0].transAxes,
                        ha="right", va="top", fontsize=8.5, family="monospace",
                        bbox=dict(fc="white", ec="0.5"))

        ref = np.array([r["raw_corr"] for r in rr])
        E = np.array([r["energy_true"] for r in rr])
        for m in ("corr", "pos", "flat"):
            y = np.array([r[m+"_corr"] for r in rr])
            g = np.isfinite(y) & np.isfinite(ref) & (ref > 0)
            axs[1].plot(E[g], 100*(y[g]/ref[g]-1), MMK[m]+"-", color=MCOL[m], ms=6, label=m)
        axs[1].axhline(0, color="k", lw=1.2)
        axs[1].set_ylabel("rispetto a raw [%]"); axs[1].grid(alpha=.3); axs[1].legend(fontsize=9)

        axs[2].plot(E, [r["pos_term"] for r in rr], "s-", color="C0", ms=6, label="POS")
        axs[2].plot(E, [r["syst_pct"] for r in rr], "d-", color="C4", ms=6,
                    label="$\\max|$trattamento $-$ raw$|$")
        axs[2].set_ylabel("[%]"); axs[2].set_xlabel("$E_{true}$ [GeV]")
        axs[2].grid(alpha=.3); axs[2].legend(fontsize=9)
        axs[2].set_xscale("log"); axs[0].set_xscale("log")
        fig.tight_layout()
        p = os.path.join(a.outdir, f"uniformita_pos_{R}ohm.png")
        fig.savefig(p, dpi=150); plt.close(fig)
        print("->", p)

    with open(os.path.join(a.outdir, "uniformita_pos_fit.csv"), "w") as fh:
        fh.write("resistance,metodo,N_MeV,err_N_MeV,S_pct,err_S,C_pct,err_C,chi2,ndf\n")
        for (R, m), (mi, n) in sorted(fitres.items()):
            fh.write(f"{R},{m},{1000*mi.values['N']:.1f},{1000*mi.errors['N']:.1f},"
                     f"{mi.values['S']:.4f},{mi.errors['S']:.4f},{mi.values['C']:.4f},"
                     f"{mi.errors['C']:.4f},{mi.fval:.2f},{n-3}\n")
    print("-> uniformita_pos_fit.csv")


if __name__ == "__main__":
    main()
