# This script will run cross-correlation between all pairs using
# xcorr.xcorr_all_ev() function.
#
# This script can be run multiple times
# A natural parallelization will be done on a 
# event basis using temproray *.lock or *.cplt files
# AND on a tag AND file size basis
# 
# OUTPUT FILE STRUCTURE:
#
#C1_04_xcorr_tag_cc_func/
#   daily  files : xcorr_tag1_tag2_000/yyyt_mm_01.h5
#                                     /yyyt_mm_02.h5
#                                     /yyyt_mm_dd.h5
#                   xcorr_tag1_tag2_001/...
#                   xcorr_tag1_tag2_XXX/...
#   master files : xcorr_tag1_tag2_XXX.h5
#
#  GROUPS inside H5 files
#   cc                       daily/hourly/... cross-correlations tables 
#   cc_nstack                order of summation for normalization 
#   in_                      input parameters passed to xcorr.xcorr_all_ev() 
#   md                       metadata for paths (lat,lon, ...)
#   md_c                     metadata for correlations
#   ref                      STACKED CROSS_CORRELATIONS
#   ref_nstack               order of summation for normalization 
#
# ...XXX.h5|- /cc  |- net1.sta1.locid |- net1.sta1.locid |- cmp1 (dataset,table)
#          |       |                  |                  |- cmp2 (dataset,table)
#          |       |                  |                  |- cmp3 (dataset,table)
#          |       |                  |                  |- ...
#          |       |                  | 
#          |       |                  |- net1.sta2.locid |- ...
#          |       |                  |- net1.stax.locid |- ...
#          |       |                  |- nety.staz.locid |- ...
#          |       | 
#          |       |- net1.sta2.locid
#          |       |- netx.stay.locid
#          |                
#          |- /ref | same as cc ... (dataset,array)
#          |- /other groups
#
########################################################
######################################################## 
########################################################
import m_pycorr.m31_xcorr_ev as xcorr 
from obspy.core import UTCDateTime
import numpy as np
in_ = {}

############### PATH(S)
in_['path']              = ['./data_4.0hz/events/IO_00','./data_4.0hz/events/IO_01','./data_4.0hz/events/IO_02','./data_4.0hz/events/IO_03','./data_4.0hz/events/IO_04'] # list of path to correlate ['../daily/tag1,', '../daily/tag2']
#in_['path']              = ['./data_4.0hz/events/glob_07'] # list of path to correlate ['../daily/tag1,', '../daily/tag2']
in_['path_out']          = './' # output path
in_['tag']               = 'IO_storm_P-PKP_PWS_2d' # output label (opt)
in_['file_size']         = 5. # [Gb] maximum size of each final h5 file, an good way to increase parallelization 

############## DATES
in_['start_time']        = -24*3600 # start correlating n sec after source time 
in_['time_win']          = 3600*2 # time window to correlate in sec
in_['time_overlap']      = 3600 # time window to correlate in sec

############## CC PARAMS
in_['cc_maxlag']         = 1800. # [s], correlations will corresponds to cc_maxlag*2*fe+1 samples 
in_['cc_cmp']            = ['ZZ','EE','NN','NE','EN']  # list of channels combination to correlate. 
in_['cc_func']           = 'ctp_coher' # "method" to compute CC :
#       'ctp_xcorr'      => Raw cross-correlation CC = tr1 . tr2*  
#       'ctp_xcorr_norm' => Energy Normalized cross-correlation CC = tr1 . tr2* / ( sqrt( sum(tr1)^2 . sum(tr1)^2 ) )
#       'ctp_coher'      => Spectral normalization,  cross-coherence CC = tr1 . tr2* / ( abs(tr1) *  abs(tr2))
in_['cc_dtype']          = 'float32' # CC precision for storage (float16,float32,float64)
in_['cc_scaling']        = 1.  # scaling factor before lowering precision 
in_['cc_tags']           = 3 # 1 : xcorr only intra-tag data,2 : xcorr only inter-tag data or 3: xcorr all  

############## STACKS
#

in_['keep_event_corr']   = True    # remove or not event corr (i.e keep only the stack)
in_['event_stack']       = False   # if more than one : stack each segments or not ?
in_['remove_event_file'] = True

## Phase weighted stacks?
in_['pws']               = True # apply pws() instead of mean()
in_['pws_timegate']      = 30. # pws params see Schimmel and Paulsen 1997
in_['pws_power']         = 4. # pws params see Schimmel and Paulsen 1997

############## MISC

in_['gzip']              = False # h5 option
in_['pp']                = [] # same as pre-processing recipe ...
#arg_filter = {
#     'type'  : 'bp' # bp, hp or lp
#    ,'f1'    : 1/50.0 # lower corner freq for bp and lp
#    ,'f2'    : 1/5.0 # higher corner freq for bp and hp
#    ,'order' : 2.0 # filter order
#    ,'taper' : 0.01 # between 0 (off) and 1 (hann window)
#}
#in_['pp_args']           = [arg_filter]

in_['use_list_xcorr'] = True  # use the list or not , default FALSE (ie compute all)
in_['list_xcorr'] = 'list_xcorr_P-PKP.txt' 

############### RUN
xcorr.xcorr_list_ev(in_)
