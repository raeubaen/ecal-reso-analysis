#!/bin/bash
# The whole chain, both selections. Reproduces run_all.sh (centroid pass) and
# run_all_hodoscope.sh (hodoscope pass) with the staged scripts.
#
#   bash pipeline/run_pipeline.sh <directory with reco_*ohm/> [output directory] [bes directory]
#
# PYTHON must be an interpreter with PyROOT (ROOT 6.40 with python 3.14 here):
#   PYTHON=/opt/homebrew/bin/python3 bash pipeline/run_pipeline.sh ...
#
# The hodoscope pass needs the curvature in crystal units of the per-run-normalised
# profiles, which is stage 3 of the centroid pass: the centroid pass runs first.
#
# Three sets of resolution points come out, the three of the slides of 10 September 2026:
#   hodoscope/10_resolution_fit.png              hodoscope cut, conservative BES subtracted
#   hodoscope/10_resolution_fit_nominal_bes.png  hodoscope cut, nominal BES subtracted
#   centroid_codiceA/10_resolution_fit.png       centroid cut of 0.182 crystals, same recipe
# plus centroid/ with the run_all.sh recipe (position correction, map systematic).
#
# The BES table of the codiceA recipe (colls_energies_summary_<R>ohm.csv) is rebuilt in
# $BES from the collimator log when COLLIMATORS points to the xlsx export and the file
# is not already there.

set -euo pipefail

BASE=${1:?directory with reco_<R>ohm/}
OUT=${2:-plot/pipeline_out}
BES=${3:-plot/bes}
PY=${PYTHON:-python3}
PY_MPL=${PYTHON_MPL:-$PY}          # an interpreter with matplotlib, for the matplotlib final plot
COLLIMATORS=${COLLIMATORS:-}
HERE=$(cd "$(dirname "$0")" && pwd)
RESISTANCES="340 400 500"

# ---------------------------------------------------------------- centroid (run_all.sh)
CEN=$OUT/centroid
mkdir -p "$CEN"
# run_all.sh drops 20592 (nominally 150 GeV, response of an 80 GeV run) in this chain
$PY "$HERE/s01_fit_dcb_per_run.py"   --base "$BASE" --outdir "$CEN" --selection centroid \
    --resistances $RESISTANCES --exclude-runs 20592
$PY "$HERE/s02_combine_runs.py"      --workdir "$CEN"
$PY "$HERE/s03_maps_and_profiles.py" --base "$BASE" --workdir "$CEN" --resistances $RESISTANCES --exclude-runs 20592
$PY "$HERE/s04_fit_parabolas.py"     --workdir "$CEN"
$PY "$HERE/s05_average_parabolas.py" --workdir "$CEN"
for MODE in run energy mean; do
  $PY "$HERE/s06_correct_amplitudes.py" --base "$BASE" --workdir "$CEN" --mode $MODE \
      --resistances $RESISTANCES --exclude-runs 20592
done
$PY "$HERE/s07_fit_corrected.py"     --workdir "$CEN"
$PY "$HERE/s08_systematics.py"       --workdir "$CEN" --besdir "$BES"
$PY "$HERE/s09_resolution_plots.py"  --workdir "$CEN"
$PY "$HERE/s10_fit_resolution.py"    --workdir "$CEN"
$PY "$HERE/s11_final_plot.py"        --workdir "$CEN"
$PY_MPL "$HERE/s11_final_plot_matplotlib.py" --workdir "$CEN" || echo "matplotlib figure skipped (no matplotlib in $PY_MPL)"

# ---------------------------------------------------------------- BES table of the codiceA recipe
if [ -n "$COLLIMATORS" ] && [ ! -f "$BES/colls_energies_summary_340ohm.csv" ]; then
  $PY "$HERE/bes_from_collimators.py" --collimators "$COLLIMATORS" \
      --timestamps "$HERE/../timestamps_runs.txt" --runs-csv "$CEN/01_dcb_per_run.csv" --outdir "$BES"
fi

# ---------------------------------------------------------------- hodoscope (run_all_hodoscope.sh)
HODO=$OUT/hodoscope
mkdir -p "$HODO"
$PY "$HERE/s01_fit_dcb_per_run.py"   --base "$BASE" --outdir "$HODO" --selection hodoscope \
    --resistances $RESISTANCES --curvature-csv "$CEN/03_profiles.csv"
$PY "$HERE/s02_combine_runs.py"      --workdir "$HODO"
$PY "$HERE/s08_systematics.py"       --workdir "$HODO" --besdir "$BES"
for BESKIND in cons nominal; do
  $PY "$HERE/s09_resolution_plots.py"  --workdir "$HODO" --bes $BESKIND
  $PY "$HERE/s10_fit_resolution.py"    --workdir "$HODO" --bes $BESKIND
  $PY "$HERE/s11_final_plot.py"        --workdir "$HODO" --bes $BESKIND
  $PY_MPL "$HERE/s11_final_plot_matplotlib.py" --workdir "$HODO" --bes $BESKIND || echo "matplotlib figure skipped"
done

# ------------------------------------------------ centroid, codice A recipe (its "cen" column)
CENA=$OUT/centroid_codiceA
mkdir -p "$CENA"
$PY "$HERE/s01_fit_dcb_per_run.py"   --base "$BASE" --outdir "$CENA" --selection centroid --recipe codiceA \
    --resistances $RESISTANCES
$PY "$HERE/s02_combine_runs.py"      --workdir "$CENA"
$PY "$HERE/s08_systematics.py"       --workdir "$CENA" --besdir "$BES"
$PY "$HERE/s09_resolution_plots.py"  --workdir "$CENA"
$PY "$HERE/s10_fit_resolution.py"    --workdir "$CENA"
$PY "$HERE/s11_final_plot.py"        --workdir "$CENA"
$PY_MPL "$HERE/s11_final_plot_matplotlib.py" --workdir "$CENA" || echo "matplotlib figure skipped"

echo "done: $CEN, $HODO and $CENA"
