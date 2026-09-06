BASE=/eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco/

# non serve, deve fare i plot il codice A
python3 drift_dcb_all.py --base $BASE --outdir plot

# non serve... deve fare i plot li fa il codice A
python3 profili_pernorm.py --base $BASE --outdir plot/profili \
    --resistances 340 400 500

# codice A
python3 resolution_hodo.py --base $BASE --outdir plot/hodo_parab --besdir plot/bes \
    --plotdir plot --resistances 340 400 500 --window parabola

# non serve... deve fare i plot li fa il codice A
python3 hodo_windows.py --base $BASE --outdir plot/hodo_parab --plotdir plot \
    --resistances 340 400 500 --window parabola
