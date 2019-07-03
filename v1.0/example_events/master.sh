#!/bin/bash

python 00_get_stations.py
python 01_config_download.py
python 02_download_data.py
python 03_preprocess_noise.py

for i in `seq 1 ${1}`;
do
        python 04_xcorr.py > stdout_04_${i}.txt &
        echo $i
        pids[${i}]=$!
done
# wait for all pids
for pid in ${pids[*]}; do
    wait $pid
done
echo "done"