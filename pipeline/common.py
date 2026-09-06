"""
What more than one stage needs: constants, event reading, the double Crystal Ball
fit, the run-combination statistics and the CSV / canvas helpers.

Everything numerical that is a FIT goes through ROOT (Minuit2 via ROOT::Fit::Fitter,
with an explicit HESSE after MIGRAD, exactly as the iminuit code did with migrad()
followed by hesse()). numpy is used only to select events and to compute plain
statistics such as means, RMS and quantiles.

Two selections share the pipeline and carry two different recipes, the ones the two
drivers of the repository use:

    centroid   -> recipe "uniforme"  (run_all.sh: uniformita_pos / uniformita_maps /
                                      resolution_final_uniforme)
    hodoscope  -> recipe "codiceA"   (run_all_hodoscope.sh: resolution_hodo.py)

The recipe decides the weights of the run combination, the systematic terms and the
N/S/C fit conventions. It is derived from the selection here, in one place, and every
stage reads it from the CSV of the previous one.
"""

import csv
import glob
import math
import os
import re
import sys

import numpy as np
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # runsets.py lives in the repository root
import runsets                                      # noqa: E402  (after the path insert)

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError      # HESSE failures are handled, not printed
ROOT.TH1.AddDirectory(False)
ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")

# ------------------------------------------------------------------ constants
SCALE = {340: 3500 / 150., 400: 1080 / 40., 500: 3340 / 100.}   # ADC per GeV, from fit.sh
FILES = {340: ("reco_340ohm", "*_merged.root"),
         400: ("reco_400ohm", "*_400_merged.root"),
         500: ("reco_500ohm", "*_500_merged.root")}
TRUE_ENERGY = {20: 20.00, 30: 30.00, 40: 39.99, 50: 49.98, 60: 59.97, 80: 79.90,
               100: 99.75, 120: 119.48, 150: 148.73, 175: 172.67, 200: 196.08,
               225: 218.82, 250: 240.76, 275: 261.77, 300: 281.74}

ETA_CENTRE, PHI_CENTRE = 18., 6.
A_TOT_MIN = 100.                 # threshold on the A_tot branch, in every selection
SYNCHROTRON_COEFF = 1.92e-7      # sigma/E [%] = coeff * E_true^2.5
CRYSTAL_WINDOW = 3               # the 3x3 matrix around (18, 6) summed for the amplitude

SELECTIONS = ("centroid", "hodoscope")
RECIPE = {"centroid": "uniforme", "hodoscope": "codiceA"}
# half-window of the position cut in crystal units: 0.2 in the uniformity chain,
# 0.182 in resolution_hodo.py (its SEL constant)
HALF_WINDOW = {"centroid": 0.2, "hodoscope": 0.182}
# runs with fewer events than this are not fitted individually
MIN_EVENTS_PER_RUN = {"codiceA": 300, "uniforme": 0}
MIN_EVENTS_POOLED = 500

TAIL_NAMES = ("alpha_l", "alpha_h", "n_l", "n_h")
DCB_PARAMETERS = TAIL_NAMES + ("mean", "sigma", "norm")
HISTOGRAM_NBINS, HISTOGRAM_LO, HISTOGRAM_HI = 8000, 0., 8000.   # fit.sh: A_tot>>h(8000,0,8000)

COLOUR = {340: ROOT.kAzure + 2, 400: ROOT.kOrange + 7, 500: ROOT.kGreen + 2}
MARKER = {340: 20, 400: 21, 500: 22}


def true_energy(nominal):
    return TRUE_ENERGY.get(int(nominal), float(nominal))


def synchrotron_pct(energy_true):
    return SYNCHROTRON_COEFF * energy_true ** 2.5


def recipe_for(selection):
    return RECIPE[selection]


# ------------------------------------------------------------------ files
def energy_of(path):
    return int(re.match(r"(\d+)", os.path.basename(path)).group(1))


def data_files(base, resistance):
    """The merged ROOT files of one resistance, sorted by nominal energy."""
    directory, pattern = FILES[resistance]
    return sorted(glob.glob(os.path.join(base, directory, pattern)), key=energy_of)


def resistance_energy_pairs(base, resistances, excluded_points=()):
    for resistance in resistances:
        for path in data_files(base, resistance):
            energy = energy_of(path)
            if (resistance, energy) in excluded_points:
                continue
            yield resistance, energy, path


def parse_excluded_points(items):
    """--exclude 340:275 -> {(340, 275)}"""
    return {tuple(int(value) for value in item.split(":")) for item in items}


# ------------------------------------------------------------------ events
ROOT.gInterpreter.Declare(r'''
// sum of the amplitudes of the crystals inside a (2*half-1)x(2*half-1) matrix
// around ieta = 18, iphi = 6: half = 2 gives the 3x3 (ieta 17-19, iphi 5-7)
double pipeline_sum_matrix(const ROOT::RVecD &amplitude,
                           const ROOT::RVec<unsigned short> &ieta,
                           const ROOT::RVec<unsigned short> &iphi, int half) {
  double total = 0.;
  for (size_t index = 0; index < amplitude.size(); ++index)
    if (std::abs((int)ieta[index] - 18) < half && std::abs((int)iphi[index] - 6) < half)
      total += amplitude[index];
  return total;
}
// first hodoscope cluster of a plane, NaN when the plane has none
double pipeline_first_or_nan(const ROOT::RVecF &positions) {
  return positions.size() ? (double)positions[0] : std::nan("");
}
''')

MATRIX_HALF = int((CRYSTAL_WINDOW + 1) / 2)


def read_events(path, amplitude="a3x3", with_hodoscope=False):
    """All the per-event quantities the pipeline uses, as numpy arrays.

    amplitude  a3x3  the 3x3 sum around (18, 6), rebuilt from the A branch, as in
                     resolution_hodo.py and drift_dcb_all.py
               atot  the A_tot branch, as in uniformita_pos.py and profili_pernorm.py
    The threshold A_tot > A_TOT_MIN is ALWAYS applied to the A_tot branch, whatever
    the amplitude fitted: that is what every flat script does.
    """
    frame = ROOT.RDataFrame("h4_reco", path)
    if amplitude == "a3x3":
        frame = frame.Define("amplitude",
                             f"pipeline_sum_matrix(A, sel_ieta, sel_iphi, {MATRIX_HALF})")
    elif amplitude == "atot":
        frame = frame.Define("amplitude", "A_tot")
    else:
        raise ValueError(f"unknown amplitude {amplitude}")
    columns = ["run", "spill", "evt", "A_tot", "amplitude", "pos_eta", "pos_phi"]
    if with_hodoscope:
        for plane in ("x1", "x2", "y1", "y2"):
            frame = frame.Define(f"hodo_{plane}", f"pipeline_first_or_nan(hodo_{plane}_pos)")
            columns += [f"hodo_{plane}", f"hodo_{plane}_nclusters"]
    arrays = frame.AsNumpy(columns)
    events = {name: np.asarray(arrays[name]) for name in columns}
    for name in ("A_tot", "amplitude", "pos_eta", "pos_phi"):
        events[name] = events[name].astype(float)
    events["run"] = events["run"].astype(np.int64)
    events["u"] = events["pos_eta"] - ETA_CENTRE
    events["v"] = events["pos_phi"] - PHI_CENTRE
    return events


def hodoscope_xy(events, yplane="y1"):
    """x = mean of the two x planes, y from the requested plane; NaN where the plane
    does not have exactly one cluster (hodoscope_calib.hodo_xy)."""
    both_x = ((events["hodo_x1_nclusters"] == 1) & (events["hodo_x2_nclusters"] == 1)
              & np.isfinite(events["hodo_x1"]) & np.isfinite(events["hodo_x2"]))
    hodo_x = np.where(both_x, 0.5 * (events["hodo_x1"] + events["hodo_x2"]), np.nan)
    one_y = (events[f"hodo_{yplane}_nclusters"] == 1) & np.isfinite(events[f"hodo_{yplane}"])
    hodo_y = np.where(one_y, events[f"hodo_{yplane}"], np.nan)
    return hodo_x, hodo_y


def runset_mask(run, dropped, kept_only):
    mask = np.ones(len(run), bool)
    if dropped:
        mask &= ~np.isin(run, dropped)
    if len(kept_only):
        mask &= np.isin(run, kept_only)
    return mask


# ------------------------------------------------------------------ double Crystal Ball
ROOT.gInterpreter.Declare(r'''
// DoubleSidedCrystalballFunction of dcb.cxx, with t in double precision as in the
// numpy version dcb_func of uniformita_pos.py / drift_dcb_all.py
double pipeline_dcb(double *x, double *par) {
  const double alpha_l = par[0], alpha_h = par[1], n_l = par[2], n_h = par[3];
  const double mean = par[4], sigma = par[5], norm = par[6];
  const double t = (x[0] - mean) / sigma;
  if (t < -alpha_l) {
    const double tail = (n_l / alpha_l) - alpha_l - t;
    return norm * std::exp(-0.5 * alpha_l * alpha_l)
                * std::pow(std::max(alpha_l / n_l * tail, 1e-12), -n_l);
  }
  if (t > alpha_h) {
    const double tail = (n_h / alpha_h) - alpha_h + t;
    return norm * std::exp(-0.5 * alpha_h * alpha_h)
                * std::pow(std::max(alpha_h / n_h * tail, 1e-12), -n_h);
  }
  return norm * std::exp(-0.5 * t * t);
}
''')

_unique_counter = [0]


def unique_name(prefix):
    _unique_counter[0] += 1
    return f"{prefix}_{_unique_counter[0]}"


def histogram_stats(counts, centres, low, high):
    """Mean and RMS of a histogram restricted to [low, high] (TH1 with SetRangeUser)."""
    inside = (centres >= low) & (centres <= high)
    selected_counts, selected_centres = counts[inside], centres[inside]
    total = selected_counts.sum()
    if total <= 0:
        return np.nan, np.nan, 0.
    mean = (selected_counts * selected_centres).sum() / total
    variance = (selected_counts * (selected_centres - mean) ** 2).sum() / total
    return mean, math.sqrt(max(variance, 0.)), total


def fit_window(values, energy, resistance):
    """The fit.sh recipe: scale*E*(0.95, 1.05), then mean +- 3 RMS twice."""
    counts, edges = np.histogram(values, bins=HISTOGRAM_NBINS, range=(HISTOGRAM_LO, HISTOGRAM_HI))
    counts = counts.astype(float)
    centres = 0.5 * (edges[:-1] + edges[1:])
    nominal = SCALE[resistance] * energy
    low, high = nominal * 0.95, nominal * 1.05
    for _ in range(2):
        mean, rms, _total = histogram_stats(counts, centres, low, high)
        if not np.isfinite(mean) or rms <= 0:
            return None
        low, high = mean - 3 * rms, mean + 3 * rms
    return low, high


def mode_window(values, energy, resistance):
    """Fallback window centred on the mode, for runs whose peak is not where the
    fit.sh recipe expects it."""
    nominal = SCALE[resistance] * energy
    near = values[(values > 0.5 * nominal) & (values < 1.3 * nominal)]
    if len(near) < 100:
        return None
    counts, edges = np.histogram(near, bins=150)
    mode = 0.5 * (edges[counts.argmax()] + edges[counts.argmax() + 1])
    core = near[np.abs(near - mode) < 0.08 * nominal]
    if len(core) < 50 or core.std() <= 0:
        return None
    return mode - 3 * core.std(), mode + 3 * core.std()


def window_is_usable(window, values):
    if window is None:
        return False
    low, high = window
    return high > low and ((values >= low) & (values <= high)).sum() >= 200


def freedman_diaconis_binwidth(values):
    if len(values) < 20:
        return 1.
    quartile_low, quartile_high = np.percentile(values, [25, 75])
    width = 2 * (quartile_high - quartile_low) / max(len(values), 1) ** (1. / 3.)
    return float(max(1., round(width)))


def fill_histogram(name, values, nbins, low, high):
    histogram = ROOT.TH1D(name, "", nbins, low, high)
    values = np.ascontiguousarray(values, dtype=float)
    if len(values):
        histogram.FillN(len(values), values, np.ones(len(values)))
    return histogram


def initial_step(value):
    """iminuit's default initial step: a tenth of the seed, 0.1 for a seed of zero."""
    return 0.1 * abs(value) if value != 0 else 0.1


def _run_fitter(function, data, seeds, limits, fixed_names):
    """MIGRAD then HESSE on a ROOT::Fit::BinData with the TF1 as model. The FitResult
    lives inside the Fitter, so its numbers are copied out before the Fitter dies."""
    fitter = ROOT.Fit.Fitter()
    wrapped = ROOT.Math.WrappedMultiTF1(function, 1)   # must outlive the fit: no copy
    fitter.SetFunction(wrapped, False)
    for index, name in enumerate(DCB_PARAMETERS):
        settings = fitter.Config().ParSettings(index)
        settings.SetValue(seeds[name])
        settings.SetStepSize(initial_step(seeds[name]))
        low, high = limits[name]
        if high is None:
            settings.SetLowerLimit(low)
        else:
            settings.SetLimits(low, high)
        if name in fixed_names:
            settings.Fix()
    fitter.Fit(data)
    hesse_ok = bool(fitter.CalculateHessErrors())
    return _copy_result(fitter.Result(), hesse_ok)


def _copy_result(result, hesse_ok):
    n_params = result.NPar()
    return dict(values=[float(result.Parameter(index)) for index in range(n_params)],
                errors=[float(result.Error(index)) for index in range(n_params)],
                chi2=float(result.Chi2()), ndf=int(result.Ndf()), valid=bool(result.IsValid()),
                hesse_ok=hesse_ok)


def fit_dcb(values, energy, resistance, fixed_tails=None):
    """Double Crystal Ball fit of one set of amplitudes: uniformita_pos.fit_dcb.

    The window comes from the fit.sh recipe (mode_window as fallback), the binning is
    Freedman-Diaconis on the events inside the window, empty bins are ignored and the
    bin error is sqrt(N) (Neyman chi2, as TH1::Fit and as the iminuit LeastSquares
    with ey = sqrt(max(y, 1)) on the non-empty bins). Three MIGRAD+HESSE rounds, each
    re-seeded from the previous one. When HESSE fails the four tail parameters are
    frozen at their fitted values and (mean, sigma, norm) are refitted: peak and sigma
    do not move, but the errors become computable.

    fixed_tails: dict of tail parameters held constant, from the pooled fit of the
    same energy and selection (the "fixed" model of --tails).

    Returns None when there is no usable window, otherwise a dict with peak, sigma,
    their errors, chi2, ndf, the tails, the window, the number of events inside it and
    the histogram + TF1 (to be written to a ROOT file or drawn).
    """
    window = fit_window(values, energy, resistance)
    if not window_is_usable(window, values):
        window = mode_window(values, energy, resistance)
    if not window_is_usable(window, values):
        return None
    low, high = window
    in_window = values[(values >= low) & (values <= high)]
    if len(in_window) < 200:
        return None
    binwidth = freedman_diaconis_binwidth(in_window)
    nbins = max(int(round((high - low) / binwidth)), 12)
    histogram = fill_histogram(unique_name("dcb_hist"), values, nbins, low, high)
    data = ROOT.Fit.BinData()               # default options: empty bins skipped
    ROOT.Fit.FillData(data, histogram)
    if data.NPoints() < 12:
        return None

    function = ROOT.TF1(unique_name("dcb_fun"), ROOT.pipeline_dcb, low, high, 7)
    for index, name in enumerate(DCB_PARAMETERS):
        function.SetParName(index, name)
    seeds = dict(alpha_l=2., alpha_h=2., n_l=2., n_h=2.,
                 mean=histogram.GetMean(), sigma=0.5 * (high - low) / 3.,
                 norm=float(histogram.GetMaximum()))
    limits = dict(alpha_l=(0.1, 10), alpha_h=(0.1, 10), n_l=(1, 10), n_h=(1, 10),
                  mean=(low, high), sigma=(0, high - low), norm=(0, None))
    fixed_names = set()
    if fixed_tails:
        seeds.update({name: fixed_tails[name] for name in TAIL_NAMES if name in fixed_tails})
        fixed_names = {name for name in TAIL_NAMES if name in fixed_tails}

    result = None
    for _ in range(3):
        result = _run_fitter(function, data, seeds, limits, fixed_names)
        seeds = dict(zip(DCB_PARAMETERS, result["values"]))
    errors_ok = result["hesse_ok"] and result["errors"][5] > 0 and result["errors"][4] > 0
    if not errors_ok:
        # tails frozen at their fitted values, only mean / sigma / norm refitted
        result_fixed = _run_fitter(function, data, seeds, limits, set(TAIL_NAMES))
        if result_fixed["hesse_ok"] and result_fixed["errors"][5] > 0 and result_fixed["errors"][4] > 0:
            result = result_fixed
    for index in range(7):
        function.SetParameter(index, result["values"][index])
        function.SetParError(index, result["errors"][index])
    tails = dict(zip(TAIL_NAMES, result["values"][:4]))
    return dict(peak=result["values"][4], err_peak=result["errors"][4],
                sigma=result["values"][5], err_sigma=result["errors"][5],
                chi2=result["chi2"], n_events=int(len(in_window)),
                n_bins=int(data.NPoints()),
                ndf=max(int(data.NPoints()) - (3 if fixed_tails else 7), 1),
                tails=tails, lo=float(low), hi=float(high),
                valid=result["valid"], histogram=histogram, function=function)


def relative_width(fit):
    """(100 * sigma/peak, its error) from the HESSE errors, as uniformita_pos.rel."""
    value = 100 * fit["sigma"] / fit["peak"]
    error = value * math.sqrt((fit["err_sigma"] / fit["sigma"]) ** 2
                              + (fit["err_peak"] / fit["peak"]) ** 2)
    return value, error


# ------------------------------------------------------------------ run combination
def weighted_mean(values, errors, weights):
    values, errors, weights = (np.asarray(item, float) for item in (values, errors, weights))
    good = np.isfinite(values) & np.isfinite(errors) & (weights > 0)
    if good.sum() == 0:
        return np.nan, np.nan
    values, errors, weights = values[good], errors[good], weights[good]
    mean = (values * weights).sum() / weights.sum()
    error = math.sqrt((weights ** 2 * errors ** 2).sum()) / weights.sum()
    return float(mean), float(error)


def weighted_scatter(values, weights):
    """Error on the weighted mean from the WEIGHTED VARIANCE of the values, which
    already contains the fit noise and the run-to-run spread. Undefined (NaN) with a
    single run."""
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    good = np.isfinite(values) & (weights > 0)
    values, weights = values[good], weights[good]
    if len(values) < 2:
        return np.nan
    sum_weights, sum_squares = weights.sum(), (weights ** 2).sum()
    effective_n = sum_weights ** 2 / sum_squares
    if effective_n <= 1:
        return np.nan
    mean = (weights * values).sum() / sum_weights
    return float(math.sqrt((weights * (values - mean) ** 2).sum() / sum_weights / (effective_n - 1)))


def syst_for_unit_chi2(values, errors):
    """Extra error s, in quadrature, such that the fit of the points to a constant has
    chi2/ndf = 1 (drift_dcb_all.syst_for_unit_chi2). Returns (s, chi2/ndf without s,
    weighted mean with s). s = 0 when the points are already compatible."""
    values, errors = np.asarray(values, float), np.asarray(errors, float)
    count = len(values)
    if count < 2 or not np.all(np.isfinite(values)) or not np.all(errors > 0):
        return np.nan, np.nan, np.nan

    def chi2_per_ndf(extra):
        weights = 1. / (errors ** 2 + extra ** 2)
        mean = (values * weights).sum() / weights.sum()
        return ((values - mean) ** 2 * weights).sum() / (count - 1)

    chi2_zero = chi2_per_ndf(0.)
    if chi2_zero <= 1:
        weights = 1. / errors ** 2
        return 0., chi2_zero, (values * weights).sum() / weights.sum()
    low, high = 0., max(values.max() - values.min(), errors.max())
    for _ in range(60):
        if chi2_per_ndf(high) <= 1:
            break
        high *= 2
    for _ in range(200):
        middle = 0.5 * (low + high)
        if chi2_per_ndf(middle) > 1:
            low = middle
        else:
            high = middle
    extra = 0.5 * (low + high)
    weights = 1. / (errors ** 2 + extra ** 2)
    return float(extra), float(chi2_zero), float((values * weights).sum() / weights.sum())


# ------------------------------------------------------------------ resolution function
def resolution_function(name, low, high):
    """sigma/E [%] = N/E (+) S/sqrt(E) (+) C with N in GeV, S and C in percent."""
    function = ROOT.TF1(name, "sqrt(pow(100*[0]/x, 2) + [1]*[1]/x + [2]*[2])", low, high)
    function.SetParNames("N", "S", "C")
    return function


def fit_graph(graph, function, fixed=(), lower_limit_zero=True):
    """TGraphErrors fit with MIGRAD + HESSE. fixed: parameter names to hold at their
    current value. Returns dict(values, errors, chi2, ndf, valid, hesse_ok) and leaves
    the fitted parameters in the TF1."""
    fitter = ROOT.Fit.Fitter()
    wrapped = ROOT.Math.WrappedMultiTF1(function, 1)   # must outlive the fit: no copy
    fitter.SetFunction(wrapped, False)
    data = ROOT.Fit.BinData()
    ROOT.Fit.FillData(data, graph)
    for index in range(function.GetNpar()):
        settings = fitter.Config().ParSettings(index)
        settings.SetValue(function.GetParameter(index))
        settings.SetStepSize(initial_step(function.GetParameter(index)))
        if lower_limit_zero:
            settings.SetLowerLimit(0.)
        if function.GetParName(index) in fixed:
            settings.Fix()
    fitter.Fit(data)
    hesse_ok = bool(fitter.CalculateHessErrors())
    result = _copy_result(fitter.Result(), hesse_ok)
    for index in range(function.GetNpar()):
        function.SetParameter(index, result["values"][index])
        function.SetParError(index, result["errors"][index])
    return result


def make_graph(x_values, y_values, y_errors=None, x_errors=None):
    x_values = np.ascontiguousarray(x_values, float)
    y_values = np.ascontiguousarray(y_values, float)
    if y_errors is None:
        y_errors = np.zeros(len(x_values))
    if x_errors is None:
        x_errors = np.zeros(len(x_values))
    return ROOT.TGraphErrors(len(x_values), x_values, y_values,
                             np.ascontiguousarray(x_errors, float),
                             np.ascontiguousarray(y_errors, float))


# ------------------------------------------------------------------ CSV
def format_value(value):
    if isinstance(value, float):
        return "nan" if not np.isfinite(value) else f"{value:.6g}"
    if isinstance(value, (np.floating,)):
        return format_value(float(value))
    if isinstance(value, (np.integer,)):
        return str(int(value))
    return str(value)


def write_csv(path, rows, columns):
    with open(path, "w") as handle:
        handle.write(",".join(columns) + "\n")
        for row in rows:
            handle.write(",".join(format_value(row.get(column, "")) for column in columns) + "\n")
    print("->", path, f"({len(rows)} rows)")


def _convert(text):
    if text == "":
        return ""
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def read_csv(path):
    """List of dicts, numbers converted; 'nan' becomes float nan."""
    with open(path) as handle:
        return [{key: _convert(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]


def require(path):
    if not os.path.exists(path):
        sys.exit(f"missing input: {path} -- run the previous stage first")
    return path


# ------------------------------------------------------------------ drawing
def style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTitleFontSize(0.045)
    ROOT.gStyle.SetLabelSize(0.04, "xyz")
    ROOT.gStyle.SetTitleSize(0.045, "xyz")
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetNumberContours(60)


def canvas_grid(name, n_panels, columns=4, panel_width=420, panel_height=340):
    rows = int(math.ceil(n_panels / columns))
    canvas = ROOT.TCanvas(name, name, columns * panel_width, max(rows, 1) * panel_height)
    canvas.Divide(columns, max(rows, 1))
    return canvas, rows


def save_canvas(canvas, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    canvas.SaveAs(path)
    print("->", path)


def text_box(lines, x_low, y_low, x_high, y_high, size=0.03, border_colour=ROOT.kGray + 2):
    box = ROOT.TPaveText(x_low, y_low, x_high, y_high, "NDC")
    box.SetFillColor(ROOT.kWhite)
    box.SetBorderSize(1)
    box.SetLineColor(border_colour)
    box.SetTextAlign(12)
    box.SetTextFont(82)
    box.SetTextSize(size)
    for line in lines:
        box.AddText(line)
    return box


def run_axis_frame(labels, y_low, y_high, title):
    """An empty histogram whose x axis carries one text label per run, to draw graphs
    of per-run quantities on top of ("P" option, x = run index)."""
    frame = ROOT.TH1F(unique_name("run_frame"), title, len(labels), -0.5, len(labels) - 0.5)
    for index, label in enumerate(labels):
        frame.GetXaxis().SetBinLabel(index + 1, str(label))
    frame.GetXaxis().LabelsOption("v")
    frame.SetMinimum(y_low)
    frame.SetMaximum(y_high)
    frame.SetStats(0)
    keep(frame)
    return frame


def keep(*objects):
    """ROOT draws pointers; Python garbage-collects them. Hold references here."""
    _kept.extend(objects)
    return objects[0] if len(objects) == 1 else objects


_kept = []
