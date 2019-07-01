import m_pycorr.m30_xcorr as xcorr 
from obspy.core import UTCDateTime

date1 = int(UTCDateTime(2018,1,1,0,0,0).timestamp)
date2 = int(UTCDateTime(2018,1,2,0,0).timestamp)

###############
###############

in_ = {}
in_['path']              = ['./data_1.0hz/daily/glob'] #,'../data_5hz/daily/set2']
in_['date']              = range(date1,date2+1,86400)
in_['hr1']               = range(0,24,4)       # if so  hour of the beginning of each segments 
in_['hr2']               = range(4,25,4)       #        hour of the end of each segment
in_['cc_maxlag']         = 1000.           # [s] 
in_['cc_cmp']            = ['ZZ']        # list des channels a correler. 
in_['cc_func']           = 'ctp_xcorr_norm' 
in_['cc_dtype']          = 'float16'    # 
in_['cc_scaling']        = 1000   # 
in_['tag']               = 'glob_conti'


in_['hr_stack']          = False   # If True, keep only one corr per day
in_['keep_daily_corr']   = True    # If True, keep daily corr in the final file 
in_['remove_daily_file'] = True   # rm or not daily files (keep them for debugging or if you plan to add dates

in_['file_size']         = 10**6           # maximum size of each h5 file ! 

in_['gzip']              = False # h5 option
in_['pws']               = True # phase weighted stack ?
in_['pws_timegate']      = 120.
in_['pws_power']         = 2.
in_['pp']                = [] # same as pre-processing recipe ...

arg_filter = {
     'type'  : 'bp' # bp, hp or lp
    ,'f1'    : 1/50.0 # lower corner freq for bp and lp
    ,'f2'    : 1/5.0 # higher corner freq for bp and hp
    ,'order' : 2.0 # filter order
    ,'taper' : 0.01 # between 0 (off) and 1 (hann window)
}

in_['pp_args']           = [arg_filter]


############### RUN
xcorr.xcorr_all(in_)

def test() : 
    pass 