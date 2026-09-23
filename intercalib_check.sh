python3 fit_dcb_per_run.py --base /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco_v3_gainSnapshotted_7x7/intercalib/1806/ --outdir out_26_intercalib/out2026_1806/ --resistances 340 --fallback-file fallback_hodo_2025.py --eta-center 18 --phi-center 6 --amplitude a3x3

python3 fit_dcb_per_run.py --base /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco_v3_gainSnapshotted_7x7/intercalib/1706/ --outdir out_26_intercalib/out2026_1706/ --resistances 340 --fallback-file fallback_hodo_2025.py --eta-center 17 --phi-center 6 --amplitude a3x3

python3 fit_dcb_per_run.py --base /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco_v3_gainSnapshotted_7x7/intercalib/1906/ --outdir out_26_intercalib/out2026_1906/ --resistances 340 --fallback-file fallback_hodo_2025.py --eta-center 19 --phi-center 6 --amplitude a3x3

python3 fit_dcb_per_run.py --base /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco_v3_gainSnapshotted_7x7/intercalib/1805/ --outdir out_26_intercalib/out2026_1805/ --resistances 340 --fallback-file fallback_hodo_2025.py --eta-center 18 --phi-center 5 --amplitude a3x3

python3 fit_dcb_per_run.py --base /eos/cms/store/group/dpg_ecal/comm_ecal/upgrade/testbeam/ECALTB_H4_Jun2026/Reco_v3_gainSnapshotted_7x7/intercalib/1807/ --outdir out_26_intercalib/out2026_1807/ --resistances 340 --fallback-file fallback_hodo_2025.py --eta-center 18 --phi-center 7 --amplitude a3x3

python3 combine_runs.py --workdir out_26_intercalib/out2026_1806/

python3 combine_runs.py --workdir out_26_intercalib/out2026_1906/

python3 combine_runs.py --workdir out_26_intercalib/out2026_1706/

python3 combine_runs.py --workdir out_26_intercalib/out2026_1805/

python3 combine_runs.py --workdir out_26_intercalib/out2026_1807/

mv out_26_intercalib/out2026_1807/02_per_run.csv out_26_intercalib/1807.csv
mv out_26_intercalib/out2026_1806/02_per_run.csv out_26_intercalib/1806.csv
mv out_26_intercalib/out2026_1805/02_per_run.csv out_26_intercalib/1805.csv
mv out_26_intercalib/out2026_1906/02_per_run.csv out_26_intercalib/1906.csv
mv out_26_intercalib/out2026_1706/02_per_run.csv out_26_intercalib/1706.csv

cd out_26_intercalib/

for f in $(ls -1 *.csv); do echo ${f::-4} $(cat $f | awk -F "," '{print $12}' | grep -v peak) $(cat $f | awk -F "," '{print $13}' | grep -v err_peak); done

cd -

