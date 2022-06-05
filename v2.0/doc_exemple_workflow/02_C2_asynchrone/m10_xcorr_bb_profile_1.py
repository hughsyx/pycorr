import m_pycorr.m30_xcorr as xcorr
from obspy.core import UTCDateTime
import glob

# compute correlations between the backbone and the profil 1 : 

date1 = int(UTCDateTime(2019,1,1,0,0,0).timestamp)
date2 = int(UTCDateTime(2019,12,31,0,0,0).timestamp)

path_list= ['data_5.0hz/daily/BB', 'data_5.0hz/daily/profil_1'] #, 'data_5.0hz/daily/IT']
in_ = {}
in_['path']              = path_list          # there should be no '/' at the end of the path
in_['date']              = [range(date1,date2+1,86400)]
in_['hr1']               = range(0,24,4)       # if so  hour of the beginning of each segments
in_['hr2']               = range(4,25,4)       #        hour of the end of each segment
in_['cc_maxlag']         = 300.               # [s]
in_['cc_cmp']            = ['ZZ']              # channel list to be correlated
in_['cc_func']           = 'ctp_coher'    #
in_['cc_dtype']          = 'float32'           # float 16 or 32
in_['cc_scaling']        = 1000                # multiply result by this constant
in_['path_out']          = ''  # prefix to output directory name
in_['out_dir']           = './'
in_['tag']               = ''              # append this to the dir name

in_['save']              = 'ref'               # [ref, day , hour]

in_['outer']  = True         # compute CC between set ?
in_['inner']  = False         # compute CC within set  ?
in_['vs_list'] = [[]]    #[range(0,2)]#,range(0,3)] [list(range(0,2)),[5]] # 2D list of virtual source indice :
in_['vs_set'] = []    #['CH']     # 1D list of corresponding net ['CH','IT']
in_['vs_inner']=False     # correlate vs with station within the same set?
in_['set_vs_all'] = []          #  list of set_set_list to be correlated against other
in_['dist']= [0, 500]        # list of min/max inter-station distance
in_['remove_daily_file'] = True  # rm or not daily files (keep them for debugging or if you plan to add dates
in_['file_size']         = 1     # maximum size of each h5 file !
in_['write_step']        = 10     # write daily correlation each N day
in_['verbose_lvl']       = 3    # set how much info are displayed. 0 = user : 1-3  = Mostly for programmer.
                                # 4 : display the full list of stations. Useful to know a station indice
in_['pws']               = False # phase weighted stack ?
in_['pws_timegate']      = 120.
in_['pws_power']         = 2.


xcorr.xcorr_all(in_)
