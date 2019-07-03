#!/bin/bash

python 00_get_stations.py
python 01_config_download.py
python 02_download_data.py
python 03_preprocess_noise.py

for i in `seq 1 2`;
do
        python 04_xcorr.py > stdout_04_${i}.txt &
        echo $i
done
