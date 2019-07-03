#!/bin/bash

python 00_get_stations.py
python 01_config_download.py

for script in 02_download_data 03_preprocess_noise 04_xcorr;
do
    for i in `seq 1 ${1}`;
    do
            python ${script}.py > log_${script}_${i}.stdout &
            echo $i
            pids[${i}]=$!
    done
    for pid in ${pids[*]}; do
        wait $pid
    done
done
python 05_plots.py
echo "done"