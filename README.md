# Energy resolution analysis — ECAL test beam H4, 2026

Energy resolution `σ/E = N/E ⊕ S/√E ⊕ C` of an ECAL 3×3 matrix exposed to electrons in
the H4 line at the SPS, for three values of the CATIA feedback resistance (340, 400 and
500 Ω, i.e. three gains). The observable is the sum of the amplitudes of the 3×3 matrix
around crystal (η 18, φ 6), fitted run by run with a double Crystal Ball; the point at
one energy is the weighted mean over the runs of σ/μ, from which the beam energy spread
(BES) and the synchrotron radiation are subtracted in quadrature.

The code to use is the **staged pipeline in `pipeline/`**: eleven scripts, one job
each, every stage saving what it computed, every fit and every plot done with ROOT
through PyROOT. It reproduces the flat scripts of the repository root as they are after
merge #2 (`resolution_hodo.py`, `run_all_hodoscope.sh`, `run_all.sh`); those are kept
as the reference they were validated against, see the last section.

This repository contains **code only**. Data, CSVs and plots stay out (`.gitignore`).

---

## 1. Requirements

* A python with **PyROOT** (ROOT ≥ 6.26) and numpy. Nothing else is needed: the events
  are read with RDataFrame, the fits use `ROOT::Fit::Fitter`, the plots are `TCanvas`.
  On this Mac that is `/opt/homebrew/bin/python3` (ROOT 6.40); on lxplus
  `source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh`.
* matplotlib, optional, for the twin of the final figure (`s11_final_plot_matplotlib.py`);
  the driver skips it when the interpreter has none.
* uproot, iminuit and matplotlib are needed only by the flat scripts, not by the pipeline.

## 2. Inputs

```
<data>/reco_340ohm/<E>_340_merged.root      (the August files are <E>_merged.root; both match)
<data>/reco_400ohm/<E>_400_merged.root
<data>/reco_500ohm/<E>_500_merged.root
<bes>/colls_energies_summary_<R>ohm.csv     BES per (resistance, energy): Energy, BES_formula, BES_cons
<bes>/rereco_<R>_withBES.csv                BES for the run_all.sh recipe only (column bes)
timestamps_runs.txt                         only to rebuild the BES table from the collimator log
```

The reference data are Ruben's merged files of 1 September 2026 on EOS,
`/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco/`
(the `BASE` of `run_all.sh`), and his BES tables at
`rgargiul.web.cern.ch/plots_ecal_mattia/bes/`. The tree is `h4_reco`; the branches used
are `run`, `spill`, `evt`, `A` with `sel_ieta`/`sel_iphi` (to rebuild the 3×3 sum),
`A_tot` (for the threshold), `pos_eta`, `pos_phi` (ECAL centroid, crystal units) and the
hodoscope `hodo_{x1,x2,y1,y2}_{nclusters,pos}`.

BES definitions (CERN-SL-Note-97-81 eq. 2, half-widths of collimators C3 and C8 in mm):
`BES_formula = √(C3² + C8²)/(27√3)` (nominal) and `BES_cons = C3/(27√3)` (lower bound,
"conservative", the one subtracted). If the table is missing,
`pipeline/bes_from_collimators.py` rebuilds it from the xlsx export of the collimator
jaws and the run time stamps; the driver does it when `COLLIMATORS` points to the xlsx.

## 3. Quick start

```bash
PYTHON=/opt/homebrew/bin/python3 bash pipeline/run_pipeline.sh <data> <out> <bes>
```

About ten minutes on a laptop, three passes, three output directories:

| directory | selection | recipe | what it is |
|---|---|---|---|
| `<out>/hodoscope/` | hodoscope window | `codiceA` | the nominal analysis (`run_all_hodoscope.sh`): slides 24 (conservative BES) and 25 (`_nominal_bes` files) of the ECAL Days talk |
| `<out>/centroid_codiceA/` | ECAL centroid, ±0.182 crystals | `codiceA` | the "cen" column of `resolution_hodo.py`: slide 26 |
| `<out>/centroid/` | ECAL centroid, ±0.2 crystals | `uniforme` | the `run_all.sh` chain with the response-surface correction and the map systematic |

The figures to look at are `11_resolution_S.png` (fit with S held at the 340 Ω value for
400 and 500 Ω, as in the slides) and `11_resolution_SC.png` (S and C held), each with a
`_mpl.png` twin drawn with matplotlib from the same CSVs. The numbers are in
`08_systematics.csv` (every term of every point), `09_resolution_points.csv` (the points
of the plot) and `10_resolution_fits.csv` (N, S, C).

Environment variables of the driver: `PYTHON` (interpreter with PyROOT), `PYTHON_MPL`
(interpreter with matplotlib, default `PYTHON`), `COLLIMATORS` (xlsx of the jaws, to
rebuild the BES table when absent).

## 4. The stages, one by one

Every stage is a script with `--help`. Stages that read events take `--base <data>`,
the others only `--workdir`, the output directory of the pass. A stage can be rerun
alone: it reads the CSVs of the previous ones and overwrites its own outputs.

| # | script | reads | writes |
|---|---|---|---|
| 1 | `s01_fit_dcb_per_run.py --base D --outdir W --selection {centroid,hodoscope}` | ROOT files | `01_dcb_per_run.csv` (peak, σ, errors, σ/μ per run, tails free and fixed, pooled fit as run 0), `01_windows.csv`, `01_dcb_fits.root`, `dcb/*.png` |
| 2 | `s02_combine_runs.py --workdir W` | 1 | `02_per_energy.csv` (weighted mean, stat, drift, vertex systematic), `02_per_run.csv`, `runs/*.png`, `drift_check_<R>ohm.png` |
| 3 | `s03_maps_and_profiles.py --base D --workdir W` | 1 + ROOT | `03_moments.csv`, `03_profiles.csv` (curvature in crystal units, the columns of `profili_pernorm.csv`), `03_parabola_centroide.csv`, `maps/`, `profiles/` |
| 4 | `s04_fit_parabolas.py --workdir W` | 3 | `04_surface_per_run.csv` |
| 5 | `s05_average_parabolas.py --workdir W` | 3 | `05_surface_mean.csv` |
| 6 | `s06_correct_amplitudes.py --base D --workdir W --mode {run,energy,mean}` | 1, 3, 4, 5 + ROOT | `06_corrected_<mode>.root`, `06_correction_<mode>.csv` |
| 7 | `s07_fit_corrected.py --workdir W` | 2, 5, 6 | `07_dcb_corrected.csv`, `07_uniformity.csv`, `uniformity_<R>ohm.png` |
| 8 | `s08_systematics.py --workdir W --besdir B` | 2, 7, BES | `08_systematics.csv` |
| 9 | `s09_resolution_plots.py --workdir W [--bes nominal]` | 8 | `09_resolution_points.csv`, `09_resolution_terms.png` |
| 10 | `s10_fit_resolution.py --workdir W [--bes nominal]` | 9 | `10_resolution_fits.csv`, `10_resolution_fits.root`, `10_resolution_fit*.png` |
| 11 | `s11_final_plot.py --workdir W [--bes nominal]`, `s11_final_plot_matplotlib.py` | 8, 9, 10 | `11_resolution_<variant>.png`, `..._mpl.png` |

Stages 3 to 7 are the position correction of the `run_all.sh` chain and run only in the
`centroid` pass. The hodoscope pass runs 1, 2, 8, 9, 10, 11 and needs
`--curvature-csv <out>/centroid/03_profiles.csv` at stage 1 (the parabola scan checks
its curvature against the one in crystal units), so the centroid pass comes first.

Options of stage 1 worth knowing:

| option | default | meaning |
|---|---|---|
| `--selection` | `centroid` | `centroid` cuts on `pos_eta`, `pos_phi`; `hodoscope` on the hodoscope window |
| `--recipe` | `auto` | `auto` = `uniforme` for the centroid, `codiceA` for the hodoscope; `--selection centroid --recipe codiceA` gives the "cen" column of `resolution_hodo.py` |
| `--amplitude` | `a3x3` | the 3×3 sum rebuilt from `A`; `atot` reads the `A_tot` branch (identical on the current files) |
| `--half` | 0.2 (`uniforme`) / 0.182 (`codiceA`) | half-window of the position cut in crystal units |
| `--tails` | `both` | DCB tails free, held at the pooled fit, or both (the "fixed" fit is written next to the free one) |
| `--runset`, `--exclude-runs`, `--exclude R:E` | `standard`, none, `340:275` | run selection, see §6 |
| `--yplane` | `y1` | hodoscope y plane |

`--bes nominal` (stages 9, 10, 11, `codiceA` only) subtracts the nominal BES instead of
the conservative one and drops the BES systematic from the error bar; the files get the
suffix `_nominal_bes`.

## 5. What the two recipes do

The recipe follows the selection and travels in the CSVs (`recipe` column), so each
stage picks the right convention by itself.

| | `uniforme` (centroid, `run_all.sh`) | `codiceA` (hodoscope and centroid 0.182, `resolution_hodo.py`) |
|---|---|---|
| weights of the run combination | number of selected events | 1/err² |
| minimum events per run | 200 inside the fit window | 300 |
| no run qualifies | point dropped | one pooled fit of all runs (`pooled = 1`) |
| BES table | `rereco_<R>_withBES.csv`; a point without BES is dropped | `colls_energies_summary_<R>ohm.csv` |
| central value | `√(σ² − BES² − sync²)` | `√(σ² − BES_cons² − sync²)` |
| error bar | stat ⊕ drift ⊕ map syst | drift ⊕ stat ⊕ vertex syst ⊕ BES syst ⊕ sync syst |
| N/S/C fit | per resistance (C held at 0.3 % at 500 Ω), then S and C common | per resistance; for 400 and 500 Ω S held at the 340 Ω value (`_S`), or S and C held (`_SC`) |

Terms, all in percent of σ/μ:

* **stat** — the larger of the error of the weighted mean and the weighted scatter of
  the per-run values (undefined with one run, then the fit error).
* **drift** — the extra error that brings the per-run σ/μ to χ²/ndf = 1 against a
  constant; zero with one run or when the runs are already compatible. It is **added to
  the error bar**, as the code does after merge #2 (the comment in the flat script says
  the opposite of what the line below it does).
* **sync** — `1.92·10⁻⁷ · E_true^2.5`, subtracted; its systematic is the shift of the
  point when it is scaled by 1.3.
* **BES** — subtracted; its systematic is the shift of the point between the
  conservative and the nominal BES.
* **vertex syst** (hodoscope) — RMS of the point over the four windows shifted by ±1 mm
  in x and in y.
* **map syst** (`uniforme`) — |σ corrected with the surface of the energy − σ corrected
  with the surface of the resistance|, propagated through the subtraction.
* **syst_tails** — |free tails − tails held at the pooled fit|, written to the CSVs for
  information and in no error bar (the flat script computes it and drops it).

Only BES and synchrotron are subtracted from the central value. Everything else goes
into the error bar. Changing this is a change of recipe, not of code.

## 6. Run selection

`runsets.py` is the single definition, imported by every stage that reads runs.

| set | runs | why |
|---|---|---|
| `FILTER_50MHZ` | 20486–20500, 20552–20599 | 50 MHz CATIA filter: a different readout condition (`--runset filter50`) |
| `HIGH_275` | 20652–20659 | the high-response population at 275 GeV |
| `OUTLIERS` | 20491 (120 GeV, +1.5 % response), 20788, 21116 (60 GeV 500 Ω, table position), 21033–21037 (all the 50 GeV runs at 500 Ω, beam spread), 21119 (the only 80 GeV run at 500 Ω, halo) | excluded by decision; the last two groups are the "excellent run" selection of the slides |

`--runset standard` (default) drops all three; the `centroid` pass also drops 20592 by
hand (nominally 150 GeV, response of an 80 GeV run), as `run_all.sh` does. The point
340 Ω 275 GeV is dropped everywhere (`--exclude 340:275`).

Runs of the standard set in the merged files (12 + 9 + 5 points). Not all of them enter
the per-run fits: a run needs 300 events inside the cut (`codiceA`) or 200 inside the
fit window (`uniforme`); the `runs` column of `02_per_energy.csv` lists the ones used.

| R | E [GeV] | runs |
|---|---|---|
| 340 | 20 | 20895–20899 |
| 340 | 30, 40, 60, 80, 100, 120, 150, 175 | one run each: 20541, 20530, 20528, 20526, 20521, 20474, 20535, 20539 |
| 340 | 200 | 20427 20428 20429 20434 (Ruben's files: 20434 only) |
| 340 | 225 | 20513 20514 20515 20517 20518 20615 20616 20617 |
| 340 | 250 | 20481 20482 20560–20566 20585 |
| 400 | 20 / 40 / 60 | 20753 / 20841–20843 / 20847–20849 |
| 400 | 80 | 20909 20911–20915 20917–20920 |
| 400 | 100 / 150 | 20769–20772 / 20780–20782 20786 20787 20789 20799–20801 |
| 400 | 200 / 225 | 20700–20702 / 20676–20681 |
| 400 | 250 | 20683 20684 20686–20689 20691–20696 20699 (Ruben's files: no 20690) |
| 500 | 30 / 40 | 21045–21047 / 21090–21099 |
| 500 | 60 / 100 / 150 | 21081 21082 / 21056–21058 / 20938 20950 20951 20953 20954 |

Eight of the twelve 340 Ω points have a single run: no drift can be estimated there.

## 7. The hodoscope window

x is the mean of the two x planes, y is plane y1 (y2 is jagged below zero); exactly one
cluster is required in each of the three planes, which keeps about 35 % of the events.
The response ⟨A_tot⟩ is profiled against each coordinate on the core of the beam and
fitted with a parabola over `[peak − h, peak + h]` for h = 5 … 10 mm; the vertex is the
median over the scan, accepted when it moves by less than 5 mm and the width by less
than 70 % across the scan. The window is vertex ± 0.182 · 22 mm. Where the scan fails,
`resolution_hodo.py` carries hand-set vertices (`BAD_PARABOLA_*` in
`pipeline/hodoscope_window.py`, 340 Ω 225 GeV in y, 400 Ω 20 GeV in y, 500 Ω 30–80 GeV
in y and 50 GeV in x); those points are drawn with an open circle. `01_windows.csv`
records the window, the vertex, the width and why the scan failed, point by point.

## 8. The flat scripts

The scripts of the repository root are the previous, monolithic version of the same
analysis and the reference the pipeline was validated against (`pipeline/README.md`,
sections on validation). They need `numpy uproot iminuit matplotlib`.

```bash
bash run_all_hodoscope.sh      # resolution_hodo.py --window parabola: the "codiceA" recipe
bash run_all.sh                # drift_dcb_all -> uniformita_pos/maps -> resolution_final_uniforme: the "uniforme" recipe
```

`resolution_hodo.py` writes `resolution_hodo.csv` (columns `hodo_*` and `cen_*`),
`systematics.csv` and `per_run.csv`, and is what produced the slides of the ECAL Days
(10 September 2026). `root/run_fit_root.sh` repeats the N/S/C fits in ROOT from the
CSVs of the `uniforme` chain.

Conventions: no single-letter variables in `pipeline/`; the repository is code only;
the working copy `energy-reso-fitter/plot/pipeline/` is kept identical to `pipeline/`
by hand (`cmp` file by file).
