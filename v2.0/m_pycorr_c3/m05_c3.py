import inspect, os, random, h5py, glob
import numpy as np
from obspy.core import UTCDateTime
import m_pycorr.m30_xcorr as xcorr 
import m_pycorr.m20_noise_processing as npr  # THIS DEPENDENCY IS TO BE REMOVED
import m_pycorr.mods.h5 as h5pc  
import m_pycorr.mods.dd as dd 
import m_pycorr.mods.lang as lang 
#import m_pycorr.mods.pkl as pkl 
import m_pycorr.mods.pyos as pyos 
#import m_pycorr.mods.s2d as s2d 
#import m_pycorr.mods.cc2da as cc2da 
from scipy.io import savemat 
import copy
#try :
#    ipdb.set_trace()
#except :
#    pass
        # TODO : 
# -date readable by matlab ? 
# -store indice of each date correlated ? => see later 
# -to use main_loop to compute dvv, .... : generate a c1 file from C3 ? 
#
#------------------------------------------------------------------------------------------
# NOTES/ TEST
# WI : both c1 do not have the same time vector ? => should be ok : except when storing the c1. we do not know to 
#      which time vector it corresponds. So we store the time vector of both h1,h2 in md_c
#
# WI sta1 is not found in h1 ? there can be more C3 than station pair metadata, 
#   in the case of one station is not find in h1 or h2. In the case the # C3 should be equal to zero?
# WI a virtual source do not exist :> their metadata is set to nan. 
# We always store the ZZ c1, no matter which compononets are c3_ed 
# 
#
#--------------------------------------------------------------------------
# TYPICAL USAGE : 
#   h1=pycorr_x.live_md.live_md('C1_01_xcorr___coher')
#   h2=pycorr_x.live_md.live_md('C1_01_xcorr___coher')
#   virtual_source  = h1.get_station_list()['id']
#   sta1 = h1.get_station_list()['id']
#   sta2 = h2.get_station_list()['id']
#   c3.init(h1,h2,sta1,sta2,vs,in_)
# 
#-----------------------------------------------------------------------------------------------
# NOTE ON WHAT WE STORE  :
# c3                         # [nwf]          : C3 of the reference C1 
# c3_mat                     # [nvs nwf]      : C3 of the reference C1  per virtual_source
# c3_daily                   # [nday nwf]     : C3 per date (c1 are stacked btw each pair of date)
# c3_daily_stacked           # [nwf]          : sum of C3 per_date 
# c3_daily_mat               # [nday nvs nwf] : C3 per date and virtual source 
# c3_daily_mat_stacked       # [nvs nwf]      : sum of C3 per date and virtual source 
#
#
#-------------------------------------------------------------------
# NOTE ON WHAT DATA/FUNCTION WE USE : 
#   h1, h2 pycorr_live_md class instance. We access : 
#     h1.md_c['tau']
#     h1.md_c['t']
#     h1.in_['cc_cmp']
#     h1.in_dir 
#     h1.read_c1_for_path()
#     h1.read_c1_for_path_between()
#
#  read_c1_for_path return a c1_trace class instance. We use the following methods :
#    cv.filter()
#    cv.get_window()
#  
#  to compute the correlation : we import m30_xcorr.py correlation function (ctp_xcorr, ctp_xcorr_norm,...)
#  to pre-process the c1 :  we import m20_noise_processing : whitening, normenv, 1bit 
#  for apodisation we use get_cos2_window
#
#
#
#-------------------------------------------------------------------------------------------
# CODE STRUCTURE : 
#....................................................................................
# init : 
#  initialize all information that determine what and how C3 will be computed. All info are stored in the dict ex:
#  In particular split couple of station to be C3ed into several sets. 
#   ex['in_']           :  user input 
#   ex['cpl']           :  list of couple of stations to be c3ed : [['CH.DAVOX', 'CH.DAVOX'],['CH.DAVOX', 'CH.DIX']]
#   ex['set']           :  list of sets = group of couple of stations to correlated. Each elment is a dict (see cset below)
#   ex['vs']            :  list of virtual source id : ['CH.DAVOX' 'CH.DIX' 'CH.EMBD' 'CH.SLE' 'CH.SULZ']
#   ex['total_size']    :  total size of the C3 that will be stored (in Go) 
#   ex['npath_per_file']:  number of path that will stored in each set*.h5 file
#
# .................................................................................
# main_loop : 
#   loop on each set of correlation and call c3_this_set : 
#     loop 1 : loop on each set of c3 to be computed in a random order so that we have 1 CPU per set 
#     loop 2 : loop on each set of c3 to be computed in a random order allowing several CPU per set
#
# .............................................................................
# c3_this_set : 
#   computed c3 for a given set of pair of station calling c3_this_path. 
#
#   0. Define the cset dict based on the info contained in ex['sets']:
#     cset['dir']    : C3__05_c3__05-10s__sw_normenv/xcorr_set_0000
#     cset['h5_pid'] : C3__05_c3__05-10s__sw_normenv/xcorr_set_0000.26263.h5
#     cset['search'] : C3__05_c3__05-10s__sw_normenv/xcorr_set_0000*.h5
#     cset['h5']     : C3__05_c3__05-10s__sw_normenv/xcorr_set_0000.h5
#     cset['cpl']    : list if station pair to be C3ed 
#
#   1. loop on all pair of station in a random order 
#      make one h5 file per pair of station using c3_this_path 
#   2. concatenate the result in a single C3_*/set_0000.h5 file
#   3. add all metadata
#   4. finally we get C3_*/set_0000.h5 
#
#.................................................................................
# c3_this_path : 
#   0. init hdf5 file "C3_*/xcorr_set_0001/group_0000/CH.DIX_CH.DIX.4032.h5"
#   1. loop on components in_['cmp'] : ZZ_ZZ, ZZ,NZ 
#     1.1 compute_c3_single_pair_of_station : get c3 for this component 
#     1.2 save the c3_ref, the c1, and c3_mat, c3_daily_* depending on in_['save_c3*']
#    2. close the file and rename it without pid.
#
#..................................................................................
# compute_c3_single_pair_of_stations(h1,h2,virtual_source,sta1,sta2,I_cmp1,I_cmp2,in_) :
#   1. loop on virtual source to compute c3_ref : 
#      read the c1 to be correlated using pycorr.live_md:h1.read_c1_between
#      compute c3 of this virtual source by calling compute_c3_for_a_vitual_source_all_type
#      get c3_mat and c3_ref by umming the contribution of all virtual source
#      
#   2. loop on date & virtual source if in_['date'] is set: 
#	   compute c3 for each date/virtual source
#	   get c3_daily, c3_daily_stacked, c3_daily_mat, c3_daily_mat_stacked 
#
#-------------------------------------------------------------------------------------------
# compute_c3_for_a_vitual_source_all_type(trace_object,maxlag_npts,in_)
#   generic function to compute C3 for a single pair of station, date, and components (ZZ,EE) but for all type (NN/PN)
#   1. prepare both trace to be correlated using their methods : 
#     1.1 : filter 
#     1.2 : preprocessing : whitening, onebit, normenv 
#     1.3 : apodisation of the beginning and end of each c1 
#   2. window both trace : set everything but selected window to 0 with apodisation 
#   3. split causal/acausal part of each c1 
#   4. loop c3 type : PP,NN,PN,NP 
#     4.1 compute xcorr 
#   return c3['PP/NN/NP/PN']

try : 
	import ipdb 
	import scipy.io 
except : 
	pass 


def init(h1,h2,sta1,sta2,virtual_source,inu) : 
	in_={}
	in_['p1']      = 5      ;         # period band used to determine
	in_['p2']      = 10     ;         # surface and coda waves window 
	in_['pp'] = {}                    # preprocessing of the c1  
	in_['pp']['normenv']  = 0     # normalize c1 by their enveloppe ?
	in_['pp']['white']    = 0         # whiten the spectrum between p1 and p2 ?
	in_['pp']['one_bit']     = 0  # 1 bit normalisation ?
	in_['norm_c3_mat_b4_defining_ref'] = True  # normalize contribution of each Vs 

	in_['win']               = 'sw' # sw, coda or sw+coda
	in_['gw']={}                    # how do we define coda window : 
	in_['gw']['coda_dp']     = 10   # coda start after end of sw + 10*p2 seconds
	in_['gw']['coda_length'] = 1000 # coda length [seconds]
	in_['gw']['sw_v1']       = 0.3  # km/s
	in_['gw']['sw_v2']       = 1    # km/s

	in_['cmp']            = ['ZZ_ZZ','ZZ_ZN'] # list of c3 compononents :'ZZ','EE'
	in_['c3_type']        = ['PP','NN']     # list of c3type : PP,NN,NP,PN
	in_['cc_func']        = 'ctp_xcorr_norm'# function used to compute the c3
	in_['maxlag']         = int(500)        # [s] should be less than c1 length
	in_['fsize']          = 1               # in Go  
	in_['dtype']          = 'float32'      
	in_['date']           = []              # which date do we correlate. [] = correlate the ref 
	in_['save_c3_mat']    = True            # keep the full c3 vs matrix ? 
	in_['save_c3_daily']  = True            # save daily reference c3 ? 
	in_['save_c3_daily_stacked']= True            # 
	in_['save_c3_daily_mat']=True           # save daily c3 vs matrix ? HUGE !
	in_['save_c3_daily_mat_stacked'] = True # save C3 daily matrix stacked 
	in_['save_c1_vs_sta'] = False                   # save all c1 + win_used. Should be FALSE !!!!
	in_['xcorr_inside']   = True            # compute C3 inside sta1 and sta2 list?
	in_['tag']            = ''              # to name the output 
	#single_source : 
	in_['dphi'] = 5 
	in_ = lang.parse_options(in_,inu)

	ex={}
	ex['dir']            = init_create_output_dir(in_)
	ex['in_']            = in_
	ex['vs']             = virtual_source 
	ex['cpl']            = init_get_list_of_station_pairs_to_be_correlated(sta1,sta2,in_['xcorr_inside']) 
	[nppf, ts]           = init_get_npath_per_h5_file(in_,len(virtual_source),h1.md_c['tau'])
	ex['npath_per_file'] = nppf 
	ex['total_size']     = ts 
	ex['set']            = init_define_sets(ex)
	init_display_welcome_message(ex,len(virtual_source))
	main_loop(ex,h1,h2)

def main_loop(ex,h1,h2) : 
	''' two loops : first put 1 process per set, then allow several process per set'''

	# randomize the set list 
	set_list = ex['set']
	random.shuffle(set_list) 

	# loop on set : iset = dict containing the prefix of the set + cpl list
	dd.dispc('1st loop','y','b')
	for iset in set_list : 
		if os.path.isdir(iset['dir']) or len(glob.glob(iset['search'])):
			continue 
		pyos.mkdir(iset['dir']) 
		c3_this_set(iset,ex['in_'],h1,h2,ex['vs'])

	# second loop on set : allow several process to work on the same set :
	dd.dispc('2nd loop','y','b')
	for iset in set_list :
		if len(glob.glob(iset['search'])) : 
			continue 
		c3_this_set(iset,ex['in_'],h1,h2,ex['vs'])

def c3_this_set(cset,in_,h1,h2,vs) :
	dd.dispc('    c3_this_set : '+cset['dir'],'c','b')
	#1. create group directories C3_*/set_0000/group_0 once and for all to limit i/o access : 
	cpl_list = cset['cpl'] 
	group_created=[] 
	for icpl in cpl_list : 
		cpl = cts_define_cpl(cset,icpl) 
		if not icpl[2] in group_created :
			pyos.mkdir(cpl['group_dir'])
			group_created.append(icpl[2])
	#2. loop on each couple of station in a random order, C3 it if it is not (being) done :
	random.shuffle(cpl_list) 
	for icpl in cpl_list : 
		cpl = cts_define_cpl(cset,icpl) 
		if os.path.isfile(cpl['h5']) or len(glob.glob(cpl['search']))>0 : 
			continue 
		c3_this_path(cpl,in_,h1,h2,vs)
	
	#3. check if we have to concatenate all station_pair h5 files => C3_*/xcorr_set_0000.h5 
	#3.1 if not all path have been C3ed : return 
	if not cts_are_all_path_done(cset) : 
		dd.dispc('      c3_this_set : not concatenating since there are still some path left','r','b')
		return 
	#3.2 if someone else is concatenating : return 
	if len(glob.glob(cset['search'])) >0 :
		dd.dispc('      c3_this_set: not concatenating since it already (being) done','r','b')
		return 
	# 3.3 : concatenate 
	cts_concatenate_this_set(cset,in_,h1,h2,vs)

def c3_this_path(cpl,in_,h1,h2,vs) : 
	#0. init file : C3_*/xcorr_set_0000/group_0000/CH.DAVOX_CH.DIX.pid.h5
	h5=h5py.File(cpl['h5_pid'],'w')
	#1. loop on components to be C3d :'ZZ_ZZ','ZZ_ZE', ... 
	for icmp in in_['cmp'] :
		icmp1 = icmp.split('_')[0]
		icmp2 = icmp.split('_')[1]
		if hasattr(h1,'in_rot') :
			I_cmp1 = np.where(h1.in_rot['cc_cmp']==str.encode(icmp1))[0][0] # peut crasher si la composant n'existe pas
			I_cmp2 = np.where(h1.in_rot['cc_cmp']==str.encode(icmp2))[0][0] # peut crasher si la composant n'existe pas
		else	:
			I_cmp1 = np.where(h1.in_['cc_cmp']==icmp1)[0][0] # peut crasher si la composant n'existe pas
			I_cmp2 = np.where(h2.in_['cc_cmp']==icmp2)[0][0] # ms c mieux comme ca 
		dset_root = '/'+cpl['sta'][0]+'/'+cpl['sta'][1]+'/'+icmp
		c3=compute_c3_single_pair_of_stations(h1,h2,vs,cpl['sta'][0],cpl['sta'][1],I_cmp1,I_cmp2,in_) 

		#1.1. write c3 of the reference 
		for ic3 in in_['c3_type'] : 
			h5.create_dataset('/c3'+dset_root+'/'+ic3, data= c3['c3'][ic3])
			#1.2 write extra data : c3_mat ,... 
			for idata in ['c3_mat', 'c3_daily','c3_daily_stacked','c3_daily_mat','c3_daily_mat_stacked'] :
				if in_['save_'+idata] : 
					h5.create_dataset('/'+idata+'/'+dset_root+'/'+ic3, data= c3[idata][ic3])

		if in_['save_c1_vs_sta'] : 
			h5.create_dataset('/c1_vs_sta'+dset_root+'/vs_sta1/trace',  data = c3['c1']['vs_sta1']['trace'])
			h5.create_dataset('/c1_vs_sta'+dset_root+'/vs_sta1/I_win_pos',data = c3['c1']['vs_sta1']['I_win_pos'])
			h5.create_dataset('/c1_vs_sta'+dset_root+'/vs_sta1/I_win_neg',data = c3['c1']['vs_sta1']['I_win_neg'])
			h5.create_dataset('/c1_vs_sta'+dset_root+'/vs_sta2/trace'    ,data = c3['c1']['vs_sta2']['trace'])		
			h5.create_dataset('/c1_vs_sta'+dset_root+'/vs_sta2/I_win_pos',data = c3['c1']['vs_sta2']['I_win_pos'])
			h5.create_dataset('/c1_vs_sta'+dset_root+'/vs_sta2/I_win_neg',data = c3['c1']['vs_sta2']['I_win_neg'])
		#1.3 write c1 if it exists :
		dset_c1 = '/c1/'+cpl['sta'][0]+'/'+cpl['sta'][1]+'/ZZ'#+icmp
		if dset_c1 not in h5 :
			tr,is_data=h1.read_c1_for_path(cpl['sta'][0],cpl['sta'][1],I_cmp=np.where(h1.in_['cc_cmp']=='ZZ')[0][0])
			if not is_data :
				tr,is_data=h2.read_c1_for_path(cpl['sta'][0],cpl['sta'][1],I_cmp=np.where(h2.in_['cc_cmp']=='ZZ')[0][0])
			if is_data :
				h5.create_dataset(dset_c1,data=tr.tr) 

	#2. close h5 file and rename it w/o pid : 
	h5.close()
	pyos.rename(cpl['h5_pid'],cpl['h5'])

def compute_c3_single_pair_of_stations(h1,h2,virtual_source,sta1,sta2,I_cmp1,I_cmp2,in_) :
	''' loop on date/ivs call c3_atomic. Handle all user options : 
		date, save_vsmat, save_daily_mat, norm_ivs b4 defining the reference 
		nwf = len(c3['time'])
		nvs = len(virtual_source)
	'''
	# c3            : ['PP'][nwf ]    : existe toujours 
	# c3_mat        : ['PP'][nvs,vwf] : existe toujours, mais pas forcement sauvegarde ? in_['save_c3_mat']
	# c3_daily      : ['PP'][nwf]     : if len(in_['date']) >0  & in_['save_c3_daily']
	# c3_daily_mat  : ['PP'][nwf]     : if in_['date'] + in_['save_c3_daily_mat']

	c3, ndate,nvs,nwf = csp_init_output(h1,virtual_source,in_)
	dd.dispc('    computing C3 for '+sta1+' - '+sta2,'c','n')
	#1. in all case : compute c3 of the reference and keep optionally the contribution of all vs :
	for kvs,ivs in enumerate(virtual_source) : 
		# 1.1 do not compute c3 if the virtual source is at one of the station :
		if ivs == sta1 or ivs ==sta2 : continue  
		# 1.2 read data : 
		tr1, is_data1 = h1.read_c1_for_path(ivs,sta1,I_cmp=I_cmp1,I_date=-1)
		tr2, is_data2 = h2.read_c1_for_path(ivs,sta2,I_cmp=I_cmp1,I_date=-1)

		#1.3 got to next virtual source if we do not have ivs->sta1 or ivs->sta2 c1
		if not is_data1 or not is_data2 : 
			continue 

		if in_['save_c1_vs_sta'] : 
			tr1_ori = tr1.tr.copy();
			tr2_ori = tr2.tr.copy();

		cv=[]
		cv.append(tr1) 
		cv.append(tr2)
		# 1.4 compute the c3 corresponding to this virtual source : 
		c3_ivs, win_used = compute_c3_for_a_vitual_source_all_type(cv,int(in_['maxlag']/h1.md_c['tau']),in_) 
		#ipdb.set_trace()
		# 1.5 store result : 
		if in_['save_c1_vs_sta'] :
			# 1st virtual source for which we have both c1 vs->sta1 and vs->sta2
			if not 'c1' in c3 : 
				c1_nwf = len(h1.md_c['t'])
				c3['c1'] = {}
				c3['c1']['vs_sta1']={}
				c3['c1']['vs_sta1']['trace'] = np.zeros((nvs,c1_nwf))
				c3['c1']['vs_sta1']['I_win_pos'] = np.zeros((nvs,2)) 
				c3['c1']['vs_sta1']['I_win_neg'] = np.zeros((nvs,2)) 
				c3['c1']['vs_sta2']={}
				c3['c1']['vs_sta2']['trace'] = np.zeros((nvs,c1_nwf))
				c3['c1']['vs_sta2']['I_win_pos'] = np.zeros((nvs,2)) 
				c3['c1']['vs_sta2']['I_win_neg'] = np.zeros((nvs,2)) 

			c3['c1']['vs_sta1']['trace'][kvs,:]    =  tr1_ori
			c3['c1']['vs_sta1']['I_win_pos'][kvs,:]= win_used[0]['pos']
			c3['c1']['vs_sta1']['I_win_neg'][kvs,:]= win_used[0]['neg']
			c3['c1']['vs_sta2']['trace'][kvs,:]    =  tr2_ori 
			c3['c1']['vs_sta2']['I_win_pos'][kvs,:]= win_used[1]['pos']
			c3['c1']['vs_sta2']['I_win_neg'][kvs,:]= win_used[1]['neg']
		for ic3 in in_['c3_type'] : 
			if in_['save_c3_mat'] :
				c3['c3_mat'][ic3][kvs,:]=c3_ivs[ic3] 
			if in_['norm_c3_mat_b4_defining_ref'] :
				if np.max(np.abs(c3_ivs[ic3])) > 0.0 :
					c3_ivs[ic3] = c3_ivs[ic3] / np.max(np.abs(c3_ivs[ic3]))
				c3_ivs[ic3][np.isnan(c3_ivs[ic3])]=0.0
			c3['c3'][ic3] += c3_ivs[ic3]
	#2. loop on date and compute c3 per date :
	if ndate > 0 :
		for kdate,idate in enumerate(in_['date'][0]) :
			date1 = in_['date'][0][kdate]
			date2 = in_['date'][1][kdate]
			#2.1 loop on virtual source : 
			for kvs,ivs in enumerate(virtual_source) : 
				c1_h1, is_data1 = h1.read_c1_for_path_between_date(ivs,sta1,date1,date2,I_cmp=I_cmp1)
				c1_h2, is_data2 = h2.read_c1_for_path_between_date(ivs,sta2,date1,date2,I_cmp=I_cmp2)
				# skip this virtual source if we did not get the c1 ivs->sta1 or ivs->sta2 
				if not is_data1 or not is_data2 : 
					continue 
				cv = [] 
				cv.append(c1_h1) 
				cv.append(c1_h2)
				c3_ivs, win_used = compute_c3_for_a_vitual_source_all_type(cv,in_['maxlag']/h1.md_c['tau'],in_) 
				# store result 
				for ic3 in in_['c3_type'] : 
					if in_['save_c3_daily_mat'] :        # [nday nvs nwf]
						c3['c3_daily_mat'][ic3][kdate, kvs,:] += c3_ivs[ic3]
					if in_['save_c3_daily_mat_stacked'] :# [nvs nwf] 
						c3['c3_daily_mat_stacked'][ic3][kvs,:] += c3_ivs[ic3]
					# normalizing before stacking on virtual sources : 
					if in_['norm_c3_mat_b4_defining_ref'] :
						c3_ivs[ic3] = c3_ivs[ic3] / np.max(np.abs(c3_ivs[ic3]))
					if in_['save_c3_daily'] :           # [nday nwf] # apres normalisation 
						c3['c3_daily'][ic3][kdate, :] += c3_ivs[ic3]
					if in_['save_c3_daily_stacked'] :   # [nwf] # APRES normalisation
						c3['c3_daily_stacked'][ic3] += c3_ivs[ic3]   
	# stocker dans le fichier de sortie les indice dans h1 et h2 pr qu'on comprenne ce qui c passe
	return c3 

def compute_c3_for_a_vitual_source_all_type(cv,maxlag_npts,in_) :
	xcorr_= getattr(xcorr,in_['cc_func'])
	win_used = []
	for icv in cv : 
		#filtering c1 and getting window: 
		icv_cpy= copy.copy(icv)
		icv_cpy.filter(in_['p1'],in_['p2'])
		icv_win = icv_cpy.get_window(in_['p2'],in_['gw'])
		if icv_win['success']==False :  # could not window sw/coda waves 
			c3={}
			for ic3 in in_['c3_type'] : # bc inter-station dist is too large :
				c3[ic3] =np.zeros((int(2*maxlag_npts+1)))
			return c3 
		#preprocessing : 
		if in_['pp']['normenv'] : 
			icv.tr=npr._normenv(icv.tr,1./icv.tau)[0]
		if in_['pp']['white'] :
			icv.tr=npr._whitening(icv.tr,1./icv.tau,{'f1':1./in_['p2'], 'f2':1./in_['p1'],'div':10})[0]
		if in_['pp']['one_bit'] : 
			icv.tr = npr._onebit(icv.tr,1./icv.tau)[0]
		#applying apodisation on both end of the signal + around the selected wave. 
		nwf= len(icv.tr)
		p2 = in_['p2']/icv.tau 
		icv.tr= icv.tr *get_cos2_window(nwf,p2*4,0,nwf,'inside')
		win_used.append({})
		if in_['win'] == 'coda'  :
			apo_pos=get_cos2_window(nwf,p2*2,icv_win['pos']['coda'][0],icv_win['pos']['coda'][1],'inside')
			apo_neg=get_cos2_window(nwf,p2*2,icv_win['neg']['coda'][0],icv_win['neg']['coda'][1],'inside')
			win_used[-1]['pos'] = [icv_win['pos']['coda'][0], icv_win['pos']['coda'][1]] 
			win_used[-1]['neg'] = [icv_win['neg']['coda'][0], icv_win['neg']['coda'][1]]
		elif in_['win'] == 'sw' :
			apo_pos = get_cos2_window(nwf,p2,icv_win['pos']['sw'][0],icv_win['pos']['sw'][1],'outside')
			apo_neg = get_cos2_window(nwf,p2,icv_win['neg']['sw'][0],icv_win['neg']['sw'][1],'outside')
			win_used[-1]['pos'] = [icv_win['pos']['sw'][0],icv_win['pos']['sw'][1]]
			win_used[-1]['neg'] = [icv_win['neg']['sw'][0],icv_win['neg']['sw'][1]]
		elif in_['win']=='sw+coda' : 
			win_used[-1]['pos'] = [icv_win['pos']['sw'][0],icv_win['pos']['coda'][1]]
			win_used[-1]['neg'] = [icv_win['neg']['sw'][0],icv_win['neg']['coda'][1]]
			apo_pos1=get_cos2_window(nwf,p2*2,icv_win['pos']['coda'][0],icv_win['pos']['coda'][1],'inside')
			apo_neg1=get_cos2_window(nwf,p2*2,icv_win['neg']['coda'][0],icv_win['neg']['coda'][1],'inside')
			apo_pos2 = get_cos2_window(nwf,p2,icv_win['pos']['sw'][0],icv_win['pos']['sw'][1],'outside')
			apo_neg2 = get_cos2_window(nwf,p2,icv_win['neg']['sw'][0],icv_win['neg']['sw'][1],'outside')
			apo_pos = apo_pos1+apo_pos2  
			apo_neg = apo_neg1+apo_neg2 
		apo = apo_pos+apo_neg 
		apo[np.where(apo>=1)]=1
		#for dbg purpose : 
		#savemat('toto_pre.mat',{'c1' :icv})
		icv.tr = icv.tr*(apo)
		#savemat('toto_post.mat',{'c1' :icv})
		#ipdb.set_trace()
	#splitting causal/acausal part and C3ing : 		
	cv0_ca = cv[0].split_ca()
	cv1_ca = cv[1].split_ca()
	c3 = {}
	for ic3 in in_['c3_type'] :
		c3[ic3] =xcorr_(cv0_ca[ic3[0]]['trace'],cv1_ca[ic3[1]]['trace'],maxlag_npts)
		c3[ic3][np.isnan(c3[ic3])]=0.0
		#trace1[np.isnan(trace1)] = 0.0
	return c3 , win_used 


#--------------------------------------------------------------------
#   
#                          CTS FUNCTIONS
#
#-----------------------------------------------------------------------
def cts_are_all_path_done(cset) : 
	''' all path are done = all h5 file exist'''
	for icpl in cset['cpl'] :
		cpl = cts_define_cpl(cset,icpl)
		if not os.path.isfile(cpl['h5']):
			return False
	return True

def cts_define_cpl(cset,icpl) : 
	pid= str(os.getpid())
	cpl = {} 
	cpl['sta']       = icpl 
	cpl['group_dir'] = cset['dir']+'/group_'+str(icpl[2]).zfill(4)
	cpl['prefix']    = cpl['group_dir']+'/'+icpl[0]+'_'+icpl[1] 
	cpl['h5']        = cpl['prefix']+'.h5'
	cpl['h5_pid']    = cpl['prefix']+'.'+pid+'.h5'
	cpl['search']    = cpl['prefix']+'.*h5'
	return cpl 


#----------------------------------------------------------------------------------------
#
#
#  compute_c3_single_pair_of_stations functions 
#
#
#-----------------------------------------------------------------------------------------
def csp_init_output(h1,virtual_source,in_) : 
	#1. get ndate,nvs,nwf :
	c3 ={}
	c3['time'] = np.arange(-in_['maxlag'],in_['maxlag']+h1.md_c['tau'],h1.md_c['tau'])

	ndate= len(in_['date'])
	if ndate == 2 :
		ndate = len(in_['date'][0])
	nwf  = len(c3['time'])
	nvs  = len(virtual_source)
	#2. initialize all output :  
	c3['c3']   ={}
	for kc3,ic3 in enumerate(in_['c3_type']) :
		#2.1 : c3['c3']['PP']   = c3 of the reference : 
		c3['c3'][ic3]=np.zeros((nwf))
		#2.2 if save_c3_mat we store the contribution of each vs :
		if in_['save_c3_mat'] :
			if kc3==0 : c3['c3_mat'] = {}
			c3['c3_mat'][ic3]=np.zeros((nvs,nwf))
		#2.3. take care of data stored if we compute c3 date per date :
		if ndate >0 :
			#2.3.1 : save 1 c3 per date 
			if in_['save_c3_daily'] : 
				if kc3 == 0 : c3['c3_daily'] = {}
				c3['c3_daily'][ic3] = np.zeros((ndate,nwf))
			#2.3.2  : save the sum of daily c3 : 
			if in_['save_c3_daily_stacked'] : 
				if kc3 == 0 : c3['c3_daily_stacked'] = {} 
				c3['c3_daily_stacked'][ic3] = np.zeros((nwf))
			#2.3.3 : save 1 c3 matrix per date with the contribution of all vs : 
			if in_['save_c3_daily_mat'] : 
				if kc3==0 : c3['c3_daily_mat'] = {}
				c3['c3_daily_mat'][ic3] = np.zeros((ndate,nvs,nwf)) 
			#2.3.3 : save 1 c3 matrix with the contribution of all vs stacked for all date :
			if in_['save_c3_daily_mat_stacked'] : 
				if kc3==0 : c3['c3_daily_mat_stacked']={}
				c3['c3_daily_mat_stacked'][ic3] = np.zeros((nvs,nwf))
	return c3,ndate,nvs,nwf

#--------------------------------------------------------------------
#
#   
#                          Correlate_This_Set FUNCTIONS
#
#
#-----------------------------------------------------------------------
def cts_are_all_path_done(cset) : 
	''' all path are done = all h5 file exist'''
	for icpl in cset['cpl'] :
		cpl = cts_define_cpl(cset,icpl)
		if not os.path.isfile(cpl['h5']):
			return False
	return True

def cts_define_cpl(cset,icpl) : 
	pid= str(os.getpid())
	cpl = {} 
	cpl['sta']       = icpl 
	cpl['group_dir'] = cset['dir']+'/group_'+str(icpl[2]).zfill(4)
	cpl['prefix']    = cpl['group_dir']+'/'+icpl[0]+'_'+icpl[1] 
	cpl['h5']        = cpl['prefix']+'.h5'
	cpl['h5_pid']    = cpl['prefix']+'.'+pid+'.h5'
	cpl['search']    = cpl['prefix']+'.*h5'
	return cpl 

def cts_concatenate_this_set(cset,in_,h1,h2,vs) : 
	# correlate_this_set : in_ :
	#                    : md  : elev, id, lat, lon 
	#                    : md_c: tau, time, time_c1  
	#                    : vs  : elev, id, lat, lon ,id 
	# first copy all c3_ref, c3_vsmat, c1 
	dd.dispc('    concatenating '+cset['h5'],'r','b')
	h5 = h5py.File(cset['h5_pid'],'w')
	h5_path = glob.glob(cset['dir']+'/*/*h5')
	for ih5 in h5_path :
		h5_in = h5py.File(ih5,'r')
		h5pc.copy_dataset_tree(h5_in,h5,'/c3')
		if 'c1_vs_sta' in h5_in : 
			h5pc.copy_dataset_tree(h5_in,h5,'/c1_vs_sta')
		if 'c1' in h5_in : 
			h5pc.copy_dataset_tree(h5_in,h5,'/c1')
		for idata in ['c3_mat','c3_daily','c3_daily_stacked','c3_daily_mat','c3_daily_mat_stacked']:
			if idata in h5_in :
				h5pc.copy_dataset_tree(h5_in,h5,'/'+idata)
		h5_in.close()

	# adding station pair and virtual sources metadata : 
	vs_dict= cts_cct_get_virtual_source_metadata(h1,h2,vs) 
	vs_dict['id'] = np.array(vs_dict['id'],dtype='S')
	md     = cts_cct_get_station_pair_metadata(h1,h2,cset) 
	md['id'] = np.array(md['id'],dtype='S')
	md_c   = cts_cct_get_md_c(h1,h2,in_)
	h5pc.copy_dict_to_group(h5,{'in_':in_})
	h5pc.copy_dict_to_group(h5,{'vs':vs_dict})
	h5pc.copy_dict_to_group(h5,{'md':md})
	h5pc.copy_dict_to_group(h5,{'md_c':md_c})
	h5pc.copy_dict_to_group(h5,{'in_h1':h1.in_})
	h5pc.copy_dict_to_group(h5,{'in_h2':h2.in_})
	h5.close()
	pyos.rename(cset['h5_pid'],cset['h5'])
	for ih5 in h5_path :
		pyos.remove(ih5)
	for igr in glob.glob(cset['dir']+'/*/') : 
		pyos.remove_dir(igr) 
	pyos.remove_dir(cset['dir'])

def cts_cct_get_md_c(h1,h2,in_) : 
	md_c = {}
	md_c['tau']   = h1.md_c['tau']
	md_c['t_c1a'] = h1.md_c['t'] 
	md_c['t_c1b'] = h2.md_c['t'] 
	md_c['t']     = np.arange(-in_['maxlag'],in_['maxlag']+md_c['tau'],md_c['tau'])
	md_c['c1_dir']= [h1.in_dir , h2.in_dir] # add this so that it is stored in the metadata
	if len(in_['date']) > 0 :
		md_c['date1'] = in_['date'][0]
		md_c['date2'] = in_['date'][1]
		#matlab compatible date : 
		md_c['date1_ordinal']=[]
		md_c['date2_ordinal']=[]
		for idate,kdate in enumerate(md_c['date1']) :
			cdate=UTCDateTime(md_c['date1'][idate]) 
			mdate1=cdate.toordinal()+(cdate.hour*3600+cdate.minute*60+cdate.second)/86400.
			cdate=UTCDateTime(md_c['date2'][idate]) 
			mdate2=cdate.toordinal()+(cdate.hour*3600+cdate.minute*60+cdate.second)/86400.
			md_c['date1_ordinal'].append(mdate1)
			md_c['date2_ordinal'].append(mdate2)
	else :
		md_c['date1']=[] 
		md_c['date2']=[]
	return md_c 

def cts_cct_get_station_pair_metadata(h1,h2,cset) : 
	''' scann each couple of station listed in cset. 
	    attempt to find their metadata in h1 and h2. if both are found : copy them. 
	    otherwise do not fill them (e.g the c3 may have been initialized, but not the path md)
	'''
	md = {} 
	md['elev']= [] 
	md['id']  = [] 
	md['lon'] = [] 
	md['lat'] = [] 
	sta_list=[] 
	sta_list.append(h1.get_station_list())
	sta_list.append(h2.get_station_list())
	for icpl in cset['cpl'] :
		I0 = np.where(icpl[0] == sta_list[0]['id'])[0]
		I1 = np.where(icpl[1] == sta_list[1]['id'])[0]
		if len(I0)>0 and len(I1)>0 :
			I0=I0[0]
			I1=I1[0]
			md['id'].append([sta_list[0]['id'][I0],sta_list[1]['id'][I1]])
			md['lat'].append([sta_list[0]['lat'][I0], sta_list[1]['lat'][I1]])
			md['lon'].append([sta_list[0]['lon'][I0], sta_list[1]['lon'][I1]])
			md['elev'].append([sta_list[0]['elev'][I0], sta_list[1]['elev'][I1]])
	return md 

def cts_cct_get_virtual_source_metadata(h1,h2,vs_id) : 
	''' loop on virtual source id. If a vs is found in both h1 and h2 => get their metadata. 
	If a virtual source has not found in either h1 or h2 => it could not be used to compute c3 
	=> put its metadata to nan and found to False'''
	vs  = {} 
	vs['id']  = []
	vs['elev']= [] 
	vs['lon'] = [] 
	vs['lat'] = [] 
	vs['found']=[] 
	sta_list=[] 
	sta_list.append(h1.get_station_list())
	sta_list.append(h2.get_station_list())
	for ivs in vs_id : 
		I0=np.where(ivs == sta_list[0]['id'])[0]	
		I1=np.where(ivs == sta_list[1]['id'])[0]
		if len(I0)==1 & len(I1)==1 :
			vs['id'].append(sta_list[0]['id'][I0])
			vs['lon'].append(sta_list[0]['lon'][I0])
			vs['lat'].append(sta_list[0]['lat'][I0])
			vs['elev'].append(sta_list[0]['elev'][I0])
			vs['found'].append(True)
		else :
			vs['id'].append(ivs)
			vs['lon'].append(np.nan)
			vs['lat'].append(np.nan)
			vs['elev'].append(np.nan)
			vs['found']=False 
	return vs 


#--------------------------------------------------------------------------------
#
#
#                          INIT FCTS 
#
#--------------------------------------------------------------------------------
def init_create_output_dir(in_) : 
	script_name = inspect.stack()[2]
	caller = script_name[1][0:-3].split('/')[-1]
	out_dir = 'C3_'+'_'+caller+'_'+in_['tag']+'_'+str(in_['p1']).zfill(2)+'-'+str(in_['p2']).zfill(2)+'s__'+in_['win']
	if in_['win']=='coda' or in_['win']=='sw+coda':
		out_dir += '_'+str(in_['gw']['coda_dp'])+'dp'
		out_dir += '-'+str(in_['gw']['coda_length'])+'s'
	if in_['pp']['normenv'] :
		out_dir+= '_normenv'
	if in_['pp']['white'] : 
		out_dir+='_white'
	if in_['pp']['one_bit'] :
		out_dir+='_1bit'
	pyos.mkdir(out_dir)
	return out_dir 
	
def init_define_sets(ex) : 
	pid=str(os.getpid())
	set=[]
	I_set = 0 
	npath_in_this_set = 0 
	#loop on all pair of stations :
	for icpl in ex['cpl'] : 
		npath_in_this_set +=1 
		if npath_in_this_set > ex['npath_per_file'] :
			npath_in_this_set =1
			I_set+=1 
		# we have to initialize a new set :
		if npath_in_this_set ==1 : #
			set.append({})
			set[-1]['dir']    = ex['dir']+'/xcorr_set_'+str(I_set).zfill(4)
			set[-1]['h5']     = set[-1]['dir']+'.h5'
			set[-1]['search'] = set[-1]['dir']+'*.h5'
			set[-1]['h5_pid'] = set[-1]['dir']+'.'+pid+'.h5'
			set[-1]['cpl']    = []
			set[-1]['gr']     = []
			ncpl_in_this_set  = 0 
			group_number      = 0 
		# adding couple of station in this set :
		set[-1]['cpl'].append([icpl[0],icpl[1],group_number])
		#set[-1]['gr'].append(group_number)
		ncpl_in_this_set += 1 
		if ncpl_in_this_set >= 1000 :
			group_number += 1 
			ncpl_in_this_set = 0 
	return set 

def init_get_npath_per_h5_file(in_,nvs,tau) :
	size_single_path=(in_['maxlag']*2.0+1.0) * len(in_['c3_type'])/tau/10.0
	if in_['save_c3_mat'] :
		size_single_path = size_single_path * nvs
	if len(in_['date'])>0 :  
		size_single_path = size_single_path * len(in_['date']) 
	if in_['dtype']=='float64'  : 
		size_single_path = size_single_path*2.0	
	if in_['dtype']=='float16'  : 
		size_single_path = size_single_path/2.0
	size_single_path=size_single_path/1000*0.038/1024 #in Go
	npath_per_file=max(1,int(in_['fsize']/size_single_path))
	return npath_per_file, size_single_path 

def init_get_list_of_station_pairs_to_be_correlated(sta1,sta2,xcorr_inside) :
	''' list all pair of stations : between sta1 and sta2. 
		If there are common stations in sta1 and sta2 we make sure that we do not output duplicate pair.
		If xcorr_inside = True, we add station pair combination within sta1 and sta2. 
		Otherwise we list only pairs with one station inside sta1 and one within sta2. 
	'''
	# list station pair with one station in sta1 and one in sta2 : 
	cpl=[]
	for ista1 in sta1 :
		for ista2 in sta2 :
			if not [ista2,ista1] in cpl and not [ista1,ista2] in cpl : 
				cpl.append([ista1,ista2])
	# list station pair with both stations in sta1 (or sta2) : 
	if xcorr_inside == True :
		for ksta1, jsta1 in enumerate(sta1) : 
			for ksta2, jsta2 in enumerate(sta1) : 
				if ksta2 < ksta1 : continue 
				if not [jsta2,jsta1] in cpl and not [jsta1,jsta2] in cpl : 
					cpl.append([jsta1,jsta2])

		for ksta1, jsta1 in enumerate(sta2) : 
			for ksta2, jsta2 in enumerate(sta2) : 
				if ksta2 < ksta1 : continue 
				if not [jsta2,jsta1] in cpl and not [jsta1,jsta2] in cpl :  
					cpl.append([jsta1,jsta2])
	return cpl 

def init_display_welcome_message(ex,nvs) : 
	ncpl  = str(len(ex['cpl']))
	nfile = str(len(ex['set']))
	ts    = str(ex['total_size'])
	fsize = str(ex['in_']['fsize'])
	dd.dispc('  computing C3 between '+ncpl+' pair of stations','y','n')
	dd.dispc('  there will be '+nfile+' h5 files of at most '+fsize+' Go = '+ts+' Go','y','n')


#------------------------------------------------------------------------
#
#
#               UTIlITY FUNCTIONS
#
#
#--------------------------------------------------------------------------
def get_cos2_window(nwf,napo,I1,I2,type='outside') :
	#voir si napo < 0 ? 
	napo = int(napo)
	I1   = int(I1)
	I2   = int(I2)
	napo = int(np.min([napo, (I2-I1)/2-1]))
	nwf  = int(nwf)
	cos2_up = np.power(np.cos(np.linspace(np.pi/2,0,napo)),2)
	cos2_dn = np.power(np.cos(np.linspace(0,np.pi/2,napo)),2)
	if type == 'inside' : 
		apo = np.zeros((nwf)) 
		try :
			gate = np.ones(int(I2-I1-napo*2))
		except : 
			ipdb.set_trace()
		apo[I1:I2]=np.concatenate((cos2_up,gate,cos2_dn))
	elif type == 'outside' :
		apo = np.zeros((int(nwf+napo+napo))) 
		gate = np.ones(int(I2-I1))
		apo[I1:I2+napo*2]=np.concatenate((cos2_up,gate,cos2_dn))
		apo = apo[napo:-napo]
	return apo
