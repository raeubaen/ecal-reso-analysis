# The staged pipeline

The analysis of the repository, split into ten scripts that do one job each and save
what they computed. It reproduces the two drivers as they are after merge #2, and the
three sets of resolution points of the ECAL Days slides of 10 September 2026:

| pass | driver / slide | selection | recipe of the systematics |
|---|---|---|---|
| `centroid`  | `run_all.sh` (drift_dcb_all → uniformita_pos / uniformita_maps → resolution_final_uniforme) | `\|pos_eta − 18\| ≤ 0.2`, `\|pos_phi − 6\| ≤ 0.2` | `uniforme` |
| `hodoscope` | `run_all_hodoscope.sh` (resolution_hodo.py, "codice A"); slide 24 with the conservative BES, slide 25 with the nominal BES (`--bes nominal` in stages 9 and 10) | vertex of the response parabola ± 0.182 · 22 mm on the hodoscope | `codiceA` |
| `centroid_codiceA` | the "cen" column of resolution_hodo.py; slide 26 | `\|pos_eta − 18\| ≤ 0.182`, `\|pos_phi − 6\| ≤ 0.182` | `codiceA` (`--selection centroid --recipe codiceA`) |

Every fit and every plot is done with ROOT through PyROOT: the double Crystal Ball,
the parabolas, the response surfaces (TMatrixD), the N/S/C fits and the canvases.
Minuit2 with an explicit HESSE after MIGRAD and iminuit's initial step sizes, so the
fits reproduce the iminuit ones of the flat scripts (same peak, sigma, errors and
chi2 to four digits on the same input, including a low-statistics run with two local
minima; see the validation section). numpy only selects events and computes plain statistics.

```bash
PYTHON=/opt/homebrew/bin/python3 bash pipeline/run_pipeline.sh <dir with reco_*ohm/> [out] [bes dir]
```

The interpreter must have PyROOT (ROOT 6.40 built against python 3.14 here). uproot,
awkward, iminuit and matplotlib are not needed. The centroid pass runs first: the
hodoscope window needs the curvature in crystal units of stage 3.

The codiceA recipe reads the BES from `colls_energies_summary_<R>ohm.csv`, a table
that is not in the repository. `bes_from_collimators.py` rebuilds it from the xlsx export
of the collimator jaws (`CMS_ECAL_Collimators_June26.xlsx`), the run time stamps
(`timestamps_runs.txt`, CERN local time) and the run list of stage 1: C3 = XCHV.022.131
and C8 = XCSV.022.386, `BES_formula = sqrt(C3² + C8²)/(27√3)`, `BES_cons = C3/(27√3)`
(slides 11–12). The driver runs it when `COLLIMATORS=<xlsx>` is set and the table is
missing; `colls_per_run.csv` lists the jaws run by run and flags the runs that sit
across a change.

`runsets.py` carries the "excellent run" selection of the slides: on top of 20491, 20788
and 21116 it drops 21033–21037 (every 50 GeV run at 500 Ω) and 21119 (the only 80 GeV run
at 500 Ω), so 500 Ω has five points as in the slides.

## The stages

| # | script | reads | writes |
|---|---|---|---|
| 1 | `s01_fit_dcb_per_run.py` | ROOT files | `01_dcb_per_run.csv` (peak, sigma, errors, sigma/mu, tails, window; per run, per window variation, tails free and fixed, plus the pooled fit as run 0), `01_windows.csv`, `01_dcb_fits.root`, `dcb/*.png` |
| 2 | `s02_combine_runs.py` | 1 | `02_per_energy.csv` (weighted mean, stat, drift, tails and vertex systematics), `02_per_run.csv` (deviations and pulls), `runs/*.png`, `drift_check_<R>ohm.png` |
| 3 | `s03_maps_and_profiles.py` | 1 + ROOT | `03_moments.csv` (M, b per run), `03_profiles.csv` (the columns of `profili_pernorm.csv`), `03_parabola_centroide.csv`, `maps/*.png`, `profiles/*.png` |
| 4 | `s04_fit_parabolas.py` | 3 | `04_surface_per_run.csv` |
| 5 | `s05_average_parabolas.py` | 3 | `05_surface_mean.csv` (one surface per energy, one per resistance) |
| 6 | `s06_correct_amplitudes.py --mode run\|energy\|mean` | 1, 3, 4, 5 + ROOT | `06_corrected_<mode>.root` (TTree of corrected amplitudes + TH1D per run), `06_correction_<mode>.csv` |
| 7 | `s07_fit_corrected.py` | 2, 5, 6 | `07_dcb_corrected.csv`, `07_uniformity.csv` (the columns of `uniformita_maps.csv`), `uniformity_<R>ohm.png` |
| 8 | `s08_systematics.py` | 2, 7, BES tables | `08_systematics.csv`, every term and its fraction of sigma/mu |
| 9 | `s09_resolution_plots.py` | 8 | `09_resolution_points.csv`, `09_resolution_terms.png` |
| 10 | `s10_fit_resolution.py` | 9 | `10_resolution_fits.csv` (codiceA: per resistance with S held at the 340 Ω value, and with S and N both held; uniforme: per resistance and common), `10_resolution_fits.root`, `10_resolution_fit*.png` |
| 11 | `s11_final_plot.py`, `s11_final_plot_matplotlib.py` | 8, 9, 10 | `11_resolution.png` and `11_resolution_mpl.png`: everything on one figure, the points, the fit curves with their boxes and the size of every term, legends away from the points. ROOT and matplotlib from the same CSVs; `_nominal_bes` twins with `--bes nominal` |

Stages 1, 3 and 6 read the events (RDataFrame, ~3 s per file). Everything else reads
CSV and takes seconds. The hodoscope pass runs 1, 2, 8, 9, 10, 11: with the hodoscope
there is no position correction. The matplotlib figure needs an interpreter with
matplotlib (`PYTHON_MPL`, default `PYTHON`); it is skipped when there is none.

`common.py` holds what more than one stage needs; `hodoscope_window.py` the response
profile, the parabola scan and the window of resolution_hodo.py.

## What the recipe decides

The recipe follows the selection (`common.RECIPE`) and travels in the CSVs, so every
stage picks the right convention by itself:

| | `uniforme` (centroid) | `codiceA` (hodoscope) |
|---|---|---|
| weights of the run combination | number of selected events | 1/err² |
| minimum events per run | none (200 inside the window) | 300 |
| pooled fit when no run qualifies | no | yes |
| BES table | `rereco_<R>_withBES.csv`, column `bes`; a point without BES is dropped | `colls_energies_summary_<R>ohm.csv`, `BES_cons` subtracted, `BES_formula` as variation |
| central value | `sqrt(σ² − BES² − sync²)`, no position correction | the same |
| error bar | stat ⊕ drift ⊕ map syst, with map syst = σ · \|s_energy − s_mean\| / σ_corr | drift ⊕ stat ⊕ vtx syst ⊕ BES syst ⊕ sync syst |
| N/S/C fit | per resistance, C held at 0.3 % at 500 Ω; then S, C common | per resistance, S held at the 340 Ω value for 400 and 500 Ω |

The drift is the same in both: the extra error that brings the per-run sigma/mu to
chi2/ndf = 1 against a constant, added to the error bar (merge #2).

## Things the flat code does that are reproduced on purpose

* `--amplitude a3x3` (default) rebuilds the 3×3 sum from the `A` branch as
  resolution_hodo.py and drift_dcb_all.py do; `atot` reads the `A_tot` branch as the
  uniformity scripts do. On the current files the two coincide exactly.
* The hodoscope window uses the vertex of the scan but a fixed width of 22 mm, and the
  hand-set vertices of `BAD_PARABOLA_*` override the scan for the points listed, as in
  resolution_hodo.py after merge #2. The crystal width W from the curvatures is still
  computed and written to `01_windows.csv`.
* `syst_tails`, the difference between free and fixed tails, is computed and written
  (resolution_hodo.py computes it and drops it) but enters no error bar.
* The ndf written for the codiceA fits is n − 3 (n − 2 at 500 Ω), as in the flat script,
  whatever is held fixed. For 400 and 500 Ω only S is held at the 340 Ω value: the "C
  (FIXED)" of the slide boxes is a label the flat script prints for every R ≠ 340, and
  the different C values in those boxes show C was free.

## Things the flat code does that are NOT reproduced

* `--window plateau` of resolution_hodo.py: after merge #2 it ends in a NameError
  (`wx` is only assigned in the parabola branch). Only the parabola window exists here.
* The `flat` method of uniformita_pos.py (reweighting to flat illumination, weighted
  DCB fits) and its binned `fit_surface`: run_all.sh only uses them for the diagnostic
  columns of `uniformita_pos.csv`. The position term `pos_term` is kept (stage 6).
* The centroid ("cen") column that resolution_hodo.py writes next to the hodoscope one:
  the centroid pass covers that selection with the run_all.sh recipe.
* Of the drift_dcb_all.py figures: the A_tot-against-event drift plot, the per-run 2D
  map grids, the A_tot-against-centroid scatter and the drift summary. Kept: the
  occupancy and <A> maps, the run centroids, the two profiles with the parabola, the
  per-run DCB panels, peak/sigma against the run.
* A point whose parabola scan fails and has no hand-set vertex made resolution_hodo.py
  crash; here it is skipped and the reason is written to `01_windows.csv`.

## Validation against the flat scripts (all three resistances, post merge #2)

* Hodoscope pass against `resolution_hodo.py --window parabola`: sigma/mu, stat, drift,
  vertex systematic, total error, number of events, number of runs and the window label
  agree to four digits at all 28 points (largest difference 0.0007 % on sigma/mu at
  400 Ω 100 GeV, one run whose free-tails fit converges slightly differently).
* Centroid pass against `uniformita_maps.py` (moments, apply, collect): raw, run,
  energy and mean corrected sigma/mu agree within 0.001 % at all points but three.
  There the double Crystal Ball with free tails has two local minima (n_l, n_h at their
  limit of 10 or not): 340 Ω 30 GeV raw (1.1641 flat, 1.1588 here), 400 Ω 250 GeV
  (0.4525 / 0.4508), 500 Ω 60 GeV energy surface (0.7184 / 0.7097, where the flat fit
  has the worse chi2/ndf, 1.65 against 1.52). On identical input arrays the two fitters
  return the same minimum.
* The per-run-normalised profiles (`03_profiles.csv`) reproduce `profili_pernorm.csv`
  to six digits.

## Against the ECAL Days slides (24, 25, 26)

Run on Ruben's reco files (the EOS directory of `run_all.sh`, merged files of
1 September 2026, done on lxplus with LCG_107) with his BES tables
(`rgargiul.web.cern.ch/plots_ecal_mattia/bes/`) and the run selection above, the three
passes give the point sets of the slides: 12, 9 and 5 points at 340, 400 and 500 Ω, the
same runs and the same event counts as his `hodo_parab/resolution_hodo.csv` at every
point. σ/E and its error agree to four digits at 23 of the 26 hodoscope points and at
24 of the 26 centroid points; the others (340 Ω 250 GeV, 400 Ω 200 GeV, 500 Ω 40 GeV, and
500 Ω 30/40 GeV for the centroid) differ by 0.003–0.005 %, one run each whose free-tails
double Crystal Ball converges to a different local minimum. The fits:

| case | R | pipeline | slide |
|---|---|---|---|
| hodoscope, conservative BES (slide 24) | 340 | N 288 ± 14, S 2.417 ± 0.479, C 0.319 ± 0.022, χ² 16.0/9 | N 289 ± 14, S 2.396 ± 0.484, C 0.320 ± 0.022, χ² 16.3/9 |
| | 400 | N 265 ± 12, C 0.347, χ² 5.5/6 | N 266 ± 12, C 0.348, χ² 5.3/6 |
| | 500 | N 253 ± 12, C 0.327, χ² 0.5/3 | N 254 ± 12, C 0.328, χ² 0.5/3 |
| centroid 0.182 (slide 26) | 340 | N 279 ± 7, S 3.315 ± 0.206, C 0.294 ± 0.016, χ² 22.4/9 | identical |
| | 400 | N 279 ± 4, C 0.330, χ² 29.6/6 | N 279 ± 4, C 0.338, χ² 29.6/6 |
| | 500 | N 268 ± 16, C 0.265, χ² 6.6/3 | N 266 ± 16, C 0.266, χ² 6.8/3 |
| hodoscope, nominal BES (slide 25) | 340 | N 287 ± 14, S 2.501 ± 0.460, C 0.310 | N 286 ± 14, S 2.574 ± 0.432, C 0.303 |

The nominal-BES points agree with Ruben's `hodo_corr_larger_bes` column to four digits,
so the residual difference of slide 25 is in how that slide's fit was made, not in the
points. On the August reco of this Mac the same chain gives points that differ by up to
0.02 % (different event counts per run, runs 20427–20429 and 20690 present).

`bes_from_collimators.py` rebuilds Ruben's table from the collimator log at every
point but three: 340 Ω 200 GeV (runs across a C8 change, 2 and 3 mm), 400 Ω 150 GeV
(log says C3 = 20 mm, table says 8) and 500 Ω 150 GeV (20 against 10). The rebuilt
table is kept in `plot/bes_reconstructed/` next to Ruben's in `plot/bes/`.

## Naming

Variables are spelled out; nothing is a single letter.
