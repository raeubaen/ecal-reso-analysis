#!/bin/bash

set -e

BASE=${1:-/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco/}
OUT=${2:-out/}
BES=${3:-bes/}
PY=${PYTHON:-python3}

mkdir -p $OUT

RESISTANCES="340 400 500"


declare -A STEPS=(
  [fit]="fit_dcb_per_run.py --base $BASE --outdir $OUT --resistances $RESISTANCES"
  [combine]="combine_runs.py --workdir $OUT"
  [syst]="systematics.py --workdir $OUT --besdir $BES"
  [reso_plots]="resolution_plots.py --workdir $OUT --bes cons"
  [fit_reso]="fit_resolution.py --workdir $OUT --bes cons"
  [final_plot]="final_plot.py --workdir $OUT --bes cons"
)

ORDER=(fit combine syst reso_plots fit_reso final_plot)

START=${1:-fit}

for step in "${ORDER[@]}"; do
  [[ "$step" == "$START" ]] && STARTED=1
  [[ $STARTED ]] && eval "$PY ${STEPS[$step]}"
done
