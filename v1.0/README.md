# PYCORR v1.0 #

Python3.7
Ambient seismic noise correlation package

# Basic philosophy
- retrieve seismic station information from FDSN catalog : http://service.iris.edu/irisws/fedcatalog/1/
- retrieve waveforms from FDSN-webservices and/or personal archive
- basic processing "on the fly": decimation, sensor response, gaps ...
- possible (pre-) processing before correlation
- xcorr with flexible parameters for optimized tomography and/or monitoring applications
- tensor rotation to retrieve RT information
- Basic toolbox to extract and plot correlations results from large output

# recommended install

*** for MacOSX ***
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh 
*** for Linux ***
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh 

bash ~/miniconda.sh
source .bashrc  
conda create -n mypy3 python=3.8
conda update -n base -c defaults conda
conda activate mypy3
conda install numpy scipy h5py 
conda install matplotlib cartopy pillow
conda install statsmodels scikit-learn
conda install obspy