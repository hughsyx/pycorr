# pycorr_dev

Development version of pycorr, an ambient seismic noise correlation package

**Basic philosophy**

Pycorr is a suite of python code to : 

1. retrieve seismic station information from FDSN catalog : http://service.iris.edu/irisws/fedcatalog/1/
2. retrieve waveforms (i.e noise records day per day, or event waveforms) from FDSN-webservices and/or personal archive. Basic processing  are applied "on the fly": decimation, sensor response, gaps ...
3. pre-process noise records : spectral whitening, combfilters, removing of earthquake
4. compute noise/coda correlations with flexible parameters that can be either optimized for tomography and/or monitoring applications
5. rotate ZNE correlations tensor rotation to retrieve RT information
6. Basic python and matlab toolbox to extract and plot correlations resulting from large output


**Recommanded installation**

To use pycorr within an anaconda environment you need   : 

  _conda create -n p37 python=3.7_

  _conda activate p37_

  _conda install ipython_ 

  _conda install numpy_ 

  _conda install scipy_

  _conda install h5py_ 

  _conda config --add channels conda-forge_

  _conda install obspy_

  _conda install ipdb_ 

  _conda install pandas_

  _conda install pytables_ 

  _conda install cartopy_

  _conda install statsmodels_

  _conda install scikit-learn_ 


**To access pycorr from any directory** 

 you can then add to your .zshrc or .bashrc file a function : 
 
`function switch_to_pycorr_2.1(){`

`conda activate p37`

`export PYTHONPATH`

`PYTHONPATH="/path_to_your_pycorr_directory/v2.1_beta/":$PYTHONPATH`

**To use the matlab functions to plot and manipulate the correlations** 
install the m_map package from : 
https://www.eoas.ubc.ca/~rich/map.html
