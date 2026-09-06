BASE=/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco/

# 1
python3 drift_dcb_all.py --base $BASE --outdir plot

mkdir plot/uniformita

# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 340 --energies 20 30 40 --exclude-runs 20592

# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 340 --energies 60 80 100 --exclude-runs 20592

# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 340 --energies 120 150 175 --exclude-runs 20592

# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 340 --energies 200 225 250 --exclude-runs 20592


# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 400 --energies 20 40 60 --exclude-runs 20592

# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 400 --energies 80 100 150 --exclude-runs 20592

# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 400 --energies 200 225 250 --exclude-runs 20592

# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 500 --energies 30 40 60 --exclude-runs 20592

# 2  (in blocks of a few energies: each call writes its own cache)
python3 uniformita_pos.py --base $BASE --besdir plot/bes/ --ooutdir plot/uniformita \
    --resistances 500 --energies 80 100 150 --exclude-runs 20592


python3 uniformita_pos.py --only-collect --resistances 340 400 500

# 3
python3 uniformita_maps.py --base $BASE --stage moments --resistances 340 400 500
python3 uniformita_maps.py --base $BASE --stage apply   --resistances 340 400 500
python3 uniformita_maps.py --base $BASE --stage collect

# 4
python3 resolution_final_uniforme.py --central raw --syst map --suffix _unif

# 5
bash root/run_fit_root.sh
