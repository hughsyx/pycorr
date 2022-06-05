# PYCORR 2.1 C3 !!
# Compute C2 between the profile 1 and 2, and save the full matrix, i.e the contribution of each virtual source. 

def keep_element_with(vs,str_) :
    indice_vs=[]
    k=0
    for ivs in vs :
        if str_ in ivs[0:3] :
            indice_vs.append(k)
        k=k+1
    vs=vs[indice_vs]    
    return vs


import m_pycorr_c3.m05_c3 as c3
from m_pycorr_c3.live.c1_md import c1_md 
from obspy.core import UTCDateTime

in_={}
in_['p1']      = 1      ;         # period band used to determine
in_['p2']      = 5      ;         # surface and coda waves window 
in_['pp'] = {}                    # preprocessing of the c1  
in_['pp']['normenv'] = 0   # normalize c1 by their enveloppe ?
in_['pp']['white']   = 0   # whiten the spectrum between p1 and p2 ?
in_['pp']['one_bit'] = 0   # 1 bit normalisation ?
in_['norm_c3_mat_b4_defining_ref'] = True  # normalize contribution of each Vs 

in_['win']               = 'sw' # sw, coda or sw+coda
in_['gw']={}                    # how do we define coda window : 
in_['gw']['coda_dp']     = 10   # coda start after end of sw + 10*p2 seconds
in_['gw']['coda_length'] = 1000 # coda length [seconds]
in_['gw']['sw_v1']       = 0.3  # km/s
in_['gw']['sw_v2']       = 1    # km/s

in_['cmp']            = ['ZZ_ZZ'] # list of c3 compononents :'ZZ','EE'
in_['c3_type']        = ['PP','NN']     # list of c3type : PP,NN,NP,PN
in_['cc_func']        = 'ctp_xcorr_norm'# function used to compute the c3
in_['maxlag']         = int(150)        # [s] should be less than c1 length
in_['fsize']          = 1               # in Go  
in_['dtype']          = 'float32'
in_['date']           = []               # which date do we correlate. [] = correlate the ref 
in_['save_c3_mat']    = True             # keep the full c3 vs matrix ? 
in_['save_c3_daily']  = False            # save daily reference c3 ? 
in_['save_c3_daily_stacked']= False      # 
in_['save_c3_daily_mat']=False           # save daily c3 vs matrix ? HUGE !
in_['save_c3_daily_mat_stacked'] = False # save C3 daily matrix stacked 
in_['save_c1_vs_sta'] = False            # save all c1 + win_used. Should be FALSE !!!!
in_['xcorr_inside']   = False            # compute C3 inside sta1 and sta2 list?
in_['tag']            = ''               # suffix appended to the output directory 


h0 = c1_md('C1_bb_profile_1');
h1 = c1_md('C1_bb_profile_2')
h2 = c1_md('C1_bb_profile_2')


# extract the list of virtual source, and the stations of the profile 1 and 2 from the C1 files : 
in_['tag'] = 'test_single_vs'
vs   = h0.get_station_list()['id']
sta1 = h1.get_station_list()['id']
sta2 = h2.get_station_list()['id']

# filter the list of station based on their name, to separate the virtual sources from the stations belonging to 
# one of the profile : 
vs = keep_element_with(vs,'10.')
sta1=keep_element_with(sta1,'01.')
sta2=keep_element_with(sta2,'02.')

# compute the C2 
c3.init(h1,h2,sta1,sta2,vs,in_)

