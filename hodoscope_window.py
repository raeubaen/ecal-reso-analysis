"""
The hodoscope window of resolution_hodo.py (post merge #2), used by stage 1 with

The response <A_tot> is profiled against the hodoscope coordinate on the core of the
beam (|A_tot / median - 1| < 0.10) and fitted with a parabola over [peak - h, peak + h]
for every h in SCAN_HALVES; the vertex is the median over the scan and the result is
accepted only if it does not depend on h (hodoscope_calib.parabola_scan, with the
post-merge thresholds).

Where the scan fails, resolution_hodo.py carries hand-set vertices (BAD_PARABOLA_*):
they override the scan for those (resistance, energy) points and the window is
labelled "<coords>-ecal_prof". A point with no parabola and no hand-set vertex made
the flat script crash; here it is skipped with the reason printed and recorded.

The parabola fits are weighted linear least squares; ROOT does them with the linear
fitter behind TGraphErrors::Fit("pol2").
"""

import numpy as np
import ROOT

import common

CORE_FRACTION = 0.10             # |A_tot / median - 1| < CORE_FRACTION defines the beam core
PROFILE_NBINS, PROFILE_MIN_PER_BIN, PROFILE_MIN_EVENTS = 40, 150, 3000
SCAN_HALVES = (5., 6., 7., 8., 9., 10.)    # fit half-widths around the maximum [mm]
MIN_FITS = 4                               # how many of them must give a maximum
VERTEX_SPREAD_MAX = 5.                     # allowed excursion of the vertex over the scan [mm]
WIDTH_SPREAD_MAX = 0.7                     # allowed relative spread of W over the scan
VERTEX_SHIFT_MM = 1.                       # +- shift of the window for the vertex systematic

COORDINATES = ("x", "y")

# hand-set vertices [mm] and half-widths [mm] of resolution_hodo.py
BAD_PARABOLA_VERTEX = {"x": {(500, 50): -6.5},
                       "y": {(400, 20): 11, (340, 225): 4, (500, 80): -6.8, (500, 30): 8,
                             (500, 40): -5, (500, 50): 2, (500, 60): -6}}
BAD_PARABOLA_HALF = {"x": {}, "y": {(400, 20): 1.98}}


def response_profile(coordinate, response):
    """<response> against the coordinate in PROFILE_NBINS bins between the 0.5th and
    the 99.5th percentile, error = RMS/sqrt(N) per bin. None when too few events."""
    finite = np.isfinite(coordinate)
    if finite.sum() < PROFILE_MIN_EVENTS:
        return None
    coordinate, response = coordinate[finite], response[finite]
    edges = np.linspace(np.percentile(coordinate, 0.5), np.percentile(coordinate, 99.5),
                        PROFILE_NBINS + 1)
    bin_index = np.clip(np.digitize(coordinate, edges) - 1, 0, PROFILE_NBINS - 1)
    centres, means, errors = [], [], []
    for index in range(PROFILE_NBINS):
        in_bin = bin_index == index
        if in_bin.sum() < PROFILE_MIN_PER_BIN:
            continue
        centres.append(coordinate[in_bin].mean())
        means.append(response[in_bin].mean())
        errors.append(response[in_bin].std(ddof=1) / np.sqrt(in_bin.sum()))
    if len(centres) < 6:
        return None
    return np.array(centres), np.array(means), np.array(errors)


def profile_peak(centres, means):
    """Position of the maximum, smoothed over three bins."""
    smoothed = np.convolve(means, np.ones(3) / 3., mode="same")
    smoothed[0], smoothed[-1] = means[0], means[-1]
    return float(centres[int(np.argmax(smoothed))])


def parabola_fit(centres, means, errors, low, high):
    """Weighted quadratic on the profile points inside [low, high]. Returns
    (vertex, relative curvature in %/mm^2, chi2/ndf, n) or None when the fit has no
    maximum inside its own range."""
    inside = (centres >= low) & (centres <= high)
    n_points = int(inside.sum())
    if n_points < 6:
        return None
    graph = common.make_graph(centres[inside], means[inside], errors[inside])
    quadratic = ROOT.TF1(common.unique_name("hodo_pol2"), "pol2", low, high)
    graph.Fit(quadratic, "QN")                 # polynomial: ROOT uses the linear fitter
    constant, linear, curvature = (quadratic.GetParameter(index) for index in range(3))
    if curvature >= 0:
        return None
    vertex = -linear / (2 * curvature)
    if not (low <= vertex <= high):
        return None
    relative_curvature = 100 * curvature / quadratic.Eval(vertex)
    return vertex, float(relative_curvature), quadratic.GetChisquare() / max(n_points - 3, 1), n_points


def parabola_scan(profile):
    """Vertex from the scan of the fit half-width. Always returns a
    dict with 'ok' and, when not ok, 'why'."""
    out = dict(ok=False, why="", vertex=np.nan, width=np.nan, vertex_spread=np.nan,
               width_spread=np.nan, chi2ndf=np.nan, n_fits=0, peak=np.nan)
    if profile is None:
        out["why"] = "no profile"
        return out
    centres, means, errors = profile
    peak = profile_peak(centres, means)
    out["peak"] = peak

    vertices, widths, chi2s = [], [], []
    for half_width in SCAN_HALVES:
        fit = parabola_fit(centres, means, errors, peak - half_width, peak + half_width)
        if fit is None:
            continue
        vertices.append(fit[0])
        widths.append(fit[1])
        chi2s.append(fit[2])
    out["n_fits"] = len(vertices)
    if len(vertices) < MIN_FITS:
        out["why"] = f"only {len(vertices)} of {len(SCAN_HALVES)} fit ranges give a maximum"
        return out
    vertices, widths = np.array(vertices), np.array(widths)
    vertex, width = float(np.median(vertices)), float(np.median(widths))
    out.update(vertex=vertex, width=width, vertex_spread=float(vertices.max() - vertices.min()),
               width_spread=float(widths.std() / width), chi2ndf=float(np.median(chi2s)))
    if out["vertex_spread"] > VERTEX_SPREAD_MAX:
        out["why"] = f"vertex moves by {out['vertex_spread']:.1f} mm across the fit ranges"
    elif out["width_spread"] > WIDTH_SPREAD_MAX:
        out["why"] = f"W varies by {100 * out['width_spread']:.0f} % across the fit ranges"
    else:
        out["ok"] = True
    return out


def hodoscope_windows(hodo_x, hodo_y, a_tot, base_mask, resistance, energy, half):
    """The (low, high) window in x and in y, with the scan diagnostics.

    Returns dict(windows={'x': (lo, hi) | None, 'y': ...}, fallback=[coords],
                 why={coord: reason}, scan={coord: scan dict}).
    """
    core = base_mask & (np.abs(a_tot / np.median(a_tot[base_mask]) - 1) < CORE_FRACTION)
    windows, why, fallback, scans = {}, {}, [], {}

    for coordinate_name in COORDINATES:
        coordinate = hodo_x if coordinate_name == "x" else hodo_y

        profile = response_profile(coordinate[core], a_tot[core])

        scan = parabola_scan(profile)

        scans[coordinate_name] = scan
        why[coordinate_name] = "" if scan["ok"] else scan["why"]
        windows[coordinate_name] = None
        if scan["ok"]:
            windows[coordinate_name] = (scan["vertex"] - half,
                                        scan["vertex"] + half)
        # the hand-set vertex overrides the scan for the points listed
        if (resistance, energy) in BAD_PARABOLA_VERTEX[coordinate_name]:
            vertex = BAD_PARABOLA_VERTEX[coordinate_name][(resistance, energy)]
            half_used = BAD_PARABOLA_HALF[coordinate_name].get((resistance, energy), half)
            windows[coordinate_name] = (vertex - half_used,
                                        vertex + half_used)
            fallback.append(coordinate_name)
    return dict(windows=windows, fallback=fallback, why=why, scan=scans)


def window_masks(hodo_x, hodo_y, base_mask, window_x, window_y):
    """Nominal cut and the four +- VERTEX_SHIFT_MM variations of resolution_hodo.py."""
    finite = base_mask & np.isfinite(hodo_x) & np.isfinite(hodo_y)

    def cut(shift_x, shift_y):
        return (finite
                & (hodo_x >= window_x[0] + shift_x) & (hodo_x <= window_x[1] + shift_x)
                & (hodo_y >= window_y[0] + shift_y) & (hodo_y <= window_y[1] + shift_y))

    shift = VERTEX_SHIFT_MM
    return {"nominal": cut(0., 0.),
            "x_low": cut(-shift, 0.), "x_high": cut(+shift, 0.),
            "y_low": cut(0., -shift), "y_high": cut(0., +shift)}


def window_label(fallback, selection_window="parabola"):
    return "+".join(fallback) + "-ecal_prof" if fallback else selection_window
