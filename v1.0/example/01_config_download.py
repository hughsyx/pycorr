# This script will create a directory containing a db.pkl file 
# This db.pkl contains all information to download the data 

import m_pycorr.m10_get_data as get_data
from obspy.core import UTCDateTime
in_                   = {}

list_year = [2018]

###############
# Standard
###############  

in_['ch']             = ['LH'] # ['ch1','ch2',...] It will try in this order. Better to start from the lower sampling rate
in_['cmp']            = ['Z'] # ['cmp1','cmp2',...]
in_['tag']            = 'glob' # subdirectory name for the output
in_['format']         = 'h5' #h5 only !
in_['station_list']   = 'stations.txt' # input station file 

###############
# if dc == 'SDS'
############### 
in_['sds_path']       = '/media/resif/validated_seismic_data/' #'/your/path/to/sds/archive'
in_['sds_metapath']   = '/media/resif/metadata/stationXML/'    #'/your/path/to/xml/metadata'
###############
# if data restriction
###############
in_['user_id']        = None # user id for fdsn client
in_['password']       = None # password for fdsn client
###############
# preprocessing
###############  
in_['pp']                = dict() # always 'dict()'
in_['pp']['tap_len']     = 15*60 # [s] dwld date1-tap_len -> date2+tap_lens but keep date1 -> date2
in_['pp']['freq']        = 1. # freq. to which data are decimated [hz] => output sampling rate
in_['pp']['glitch']      = False # remove glitch ?  
in_['pp']['f_prefilt']   = (0.005, in_['pp']['freq'] / 2 - 0.05 * in_['pp']['freq']) # prefilt before sensor reponse deconvolution, and after decimation
in_['pp']['remove_resp'] = True # True or False. If False, it will try to download resp from clients but DO NOT try to read from SDS archive
    
###############
# QC
###############  
in_['qc']                  = dict() # always 'dict()'
in_['qc']['max_nsegments'] = 12 # trash trace if more than x segments
in_['qc']['min_duration']  = 1 # trash trace if length than x sec
in_['qc']['timeout']       = 300 # walltime in sec for processing 1 trace

############### RUN
for iyear in list_year : # loop over a range of dates
    in_['day1']=UTCDateTime(iyear,1,1)
    in_['day2']=UTCDateTime(iyear,1,2)
    get_data.configure(in_)
