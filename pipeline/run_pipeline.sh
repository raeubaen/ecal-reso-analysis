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

set -euo pipefail

BASE=${1:?directory with reco_<R>ohm/}
OUT=${2:-plot/pipeline_out}
BES=${3:-plot/bes}
PY=${PYTHON:-python3}
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

# ---------------------------------------------------------------- hodoscope (run_all_hodoscope.sh)
HODO=$OUT/hodoscope
mkdir -p "$HODO"
$PY "$HERE/s01_fit_dcb_per_run.py"   --base "$BASE" --outdir "$HODO" --selection hodoscope \
    --resistances $RESISTANCES --curvature-csv "$CEN/03_profiles.csv"
$PY "$HERE/s02_combine_runs.py"      --workdir "$HODO"
$PY "$HERE/s08_systematics.py"       --workdir "$HODO" --besdir "$BES"
$PY "$HERE/s09_resolution_plots.py"  --workdir "$HODO"
$PY "$HERE/s10_fit_resolution.py"    --workdir "$HODO"

echo "done: $CEN and $HODO"
