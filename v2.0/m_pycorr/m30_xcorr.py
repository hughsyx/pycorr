################################################
# m30_xcorr.py 
# Laurent Stehly (UGA)
# Pierre Boue (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################

import h5py, glob, os , inspect, pickle , random , ipdb , datetime
import numpy as np, pandas as pd
from numpy import array as array 
from scipy import fftpack, io as io, signal as signal
from obspy.core import UTCDateTime

from m_pycorr.mods import h5 ,lang,dd,pyos 

import m_pycorr.m30_xcorr_scheduler as sch 

# TODO :
# - scaling 
# - separe ce qui est commun a toutes les fichiers de ce qui est specifique ?
# -- on pourrait instancier un cts_main avec tout ce qui est commun et juste ensuite copier les dictionnaires dans le 
#currant ?
# -- write step : ecriture horaire ou journaliere
# -- verifier normalization des stack : en [S] : compter en H, diviser a la fin. 
# -- changer hour en win dans le fichier de config. 
#
# RMQ : 
# nstack = nombre de jour (et non pas en [seconde])
# scaling non implemente 
# write_step = on calcule les cc jours par jour, et on les ecrit tous les N jours
def xcorr_all(inu) :
	date1 = int(UTCDateTime(2016,1,1,0,0,0).timestamp)
	date2 = int(UTCDateTime(2016,1,5,0,0).timestamp)

	in_ = {}
	in_['out_dir']='./'
	in_['path']              = ['data_5.0hz/daily/CH', 'data_5.0hz/daily/IT'] # there should be no '/' at the end 
	in_['date']              = [range(date1,date2+1,86400)] 
	in_['hr1']               = range(0,24,4)       # if so  hour of the beginning of each segments 
	in_['hr2']               = range(4,25,4)       #        hour of the end of each segment
	in_['cc_maxlag']         = 3600.               # [s] 
	in_['cc_cmp']            = ['ZZ']              # channel list to be correlated 
	in_['cc_func']           = 'ctp_xcorr_norm'    #
	in_['cc_dtype']          = 'float32'           # float 16 or 32 
	in_['cc_scaling']        = 1000                # multiply result by this constant 
	in_['path_out']          = ''                  # prefix to output directory name
	in_['tag']               = 'test'              # append this to the dir name
	in_['save']              = 'ref'               # [ref, day , hour]

	in_['outer']  = True          # compute CC between set ? 
	in_['inner']  = True          # compute CC within set  ?
	in_['vs_list'] = [] #[range(0,2)]#,range(0,3)] [list(range(0,2)),[5]] # 2D list of virtual source indice :
	in_['vs_set'] = []  #['CH']     # 1D list of corresponding net ['CH','IT']
	in_['vs_inner'] = False 
	in_['set_vs_all'] = []          #  list of set_set_list to be correlated against other 
	in_['dist']= [0, 1000000]        # list of min/max inter-station distance 
	in_['lat']=[-90,90] 
	in_['lon']=[-180,180]
	in_['remove_daily_file'] = True  # rm or not daily files (keep them for debugging or if you plan to add dates
	in_['file_size']         = 1     # maximum size of each h5 file ! 
	in_['write_step']        = 2     # write daily correlation each N day
	in_['verbose_lvl']       = 3     # set how much info are displayed. 0 = user : 1-3  = Mostly for programmer. 

	#in_['pws']               = False # not implemented  ?
	#in_['pws_timegate']      = 120.
	#in_['pws_power']         = 2.
	in_ = lang.parse_options(in_,inu)

	out_dir = get_and_create_directory_name(in_) 
	# call the scheduler to build db_sta.h5, db_cc.h5 and all db_cc_...file 
	sch.schedule_all(in_,out_dir)

	# get the list of station from db_sta.h5
	db_sta = DB_sta(out_dir)
	#get the list of (shuffled) file indice to be computed : 
	I_file = get_file_indice(out_dir)
	for ifile in I_file : 
		# set cts.I_set1,I_set2 and cts.name which contains the name of output dir 
		cts = CTS(out_dir,ifile)
		if cts.is_already_done() : 
			continue 
		cts.create_output_dir_and_lock_file()
		dd.dispc('    '+str(datetime.datetime.now())[0:-7]+' :    lets go !','c','n')
		# set cts.ex and initialize cts.h5 matrix for ref, ref_nstack and optionnaly cc_nstack
		cts.init(in_,db_sta.fe)
		# intialize the tmp h5 file with md, md_c and in_ and optionaly the tree of daily CC with resizable dataset
		cts.init_h5_file_and_get_mdc_I1(db_sta,in_)
		# loop on group of days [0-2, 2-4, 4-6,...]. At the end of each group of days we write the daily CC 
		for iday_group in range(0,len(cts.ex['I_day'][0,:])) :
			# init the matrix of daily CC for this group of days :
			if in_['save']== 'hour' or in_['save'] == 'day' :
				ntr = cts.ex['I_day'][1,iday_group] - cts.ex['I_day'][0,iday_group]
				if in_['save']=='hour' : 
					ntr = ntr * cts.ex['nhr']
				cts.h5['cc']=np.zeros((cts.ex['ncmp'],cts.ex['npath'],ntr,cts.ex['nwf']))
			# loop on days of this group of days, kdays is used to write into cts.h5['cc'] : 
			kday = -1 
			for iday in range (cts.ex['I_day'][0,iday_group],cts.ex['I_day'][1,iday_group]) :
				#print(iday) 
				kday=kday+1 
				# attempt to open data_1.0hz/daily/CH/2016/day_001.h5 files
				h5a, h5b = cts.open_noise_file(iday,db_sta)
				if h5a == False : continue
				if h5b == False : continue
				# We found both files so we loop on cmp ZZ,ZN and on station pairs 
				for icmp in range(0,len(in_['cc_cmp'])) :
					for icpl in range(0,cts.ex['npath']) :
						trace0, trace1 = cts.read_data(h5a,h5b,in_['cc_cmp'][icmp],icpl,db_sta) 
						if len(trace0) == 0 : continue 
						if len(trace1) == 0 : continue 
						# compute correlation 
						maxlag  = int(in_['cc_maxlag']*db_sta.fe)
						cc_func = in_['cc_func']
						dtype   = in_['cc_dtype']
						# AJOUTER SCALING 
						cc_hr, ncorr=correlate_this_path(trace0,trace1,maxlag,cts.I1,cts.I2,cc_func,dtype)
						# update the reference correlations, nstack
						cc_hrsum = cc_hr.sum(axis=0)
						if len(cc_hrsum[np.isnan(cc_hrsum)]) > 0  :
							ipdb.set_trace()

						cts.h5['ref'][icmp,icpl,:] += cc_hr.sum(axis=0)
						cts.h5['ref_nstack'][icmp,icpl] += 1 
						# optionnaly update cc_nstack and cc :
						if in_['save']=='day':
							cts.h5['cc'][icmp,icpl,kday,:]=cc_hr.sum(axis=0)
							cts.h5['cc_nstack'][icmp,icpl,iday]+=1	
						if in_['save']=='hour' : #kday = numero de jour au sein de ce group de jour
								Itr1 =  cts.ex['nhr']*kday
								Itr2 =  Itr1 + cts.ex['nhr']
								Itr3 =  cts.ex['nhr']*iday
								Itr4 =  Itr3 + cts.ex['nhr']-1
								cts.h5['cc'][icmp,icpl,Itr1:Itr2,:]=cc_hr#.sum(axis=0)
								cts.h5['cc_nstack'][icmp,icpl,Itr3:Itr4]=1						

				h5a.close()
				h5b.close()
			if in_['save']== 'day' :
				I1 = cts.ex['I_day'][0,iday_group]
				I2 = cts.ex['I_day'][1,iday_group]
			if in_['save']=='hour' :
				I1 = cts.ex['I_day'][0,iday_group]*cts.ex['nhr']
				I2 = cts.ex['I_day'][1,iday_group]*cts.ex['nhr']
			if in_['save'] == 'hour' or in_['save']=='day' :
				cts.write_cc(db_sta,in_['cc_cmp'],I1,I2) # on laisse cts.I1, vue que ca va change plus tard

		# end of day group everything is finished :-) 
		cts.write_all(db_sta,in_['cc_cmp'])
		cts.clean_all()



###########################################################################################
##########          ####    ############            #####            #####            ####
##########     #########    ############    ####    #####    #############    ############
#########     ##########    ############    ####    #####    #############    ############
#########     ##########    ############            #####            #####            ####
#########      #########    ############    ####    #############    #############    ####
##########          ####    ############    ####    #############    #############    ####
###########         ####            ####    ####    #####            #####            ####
##########################################################################################
##########################################################################################

class CTS :
	def __init__(self,out_dir,ifile) :
		# get the set indice from db_cc : 
		set_ = pd.read_hdf(out_dir+'/db_cc.h5','/set',start=ifile,stop=ifile+1)
		self.I_set1 = set_['I_set1'].iat[0]
		self.I_set2 = set_['I_set2'].iat[0]

		#get filenames for this set of cc :		
		pid  = str(os.getpid())
		root = set_['name'].iat[0] + '/'+ set_['name'].iat[0].split('/')[-1]
		self.name = {}
		self.name['out_dir'] = out_dir 
		self.name['dir']    = out_dir+'/'+set_['name'].iat[0]
		self.name['lock']   = out_dir+'/'+root+'.'+pid+'.lock'
		self.name['h5_tmp'] = out_dir+'/'+root+'.'+pid+'.h5'
		self.name['search'] = out_dir+'/'+root+'.*'
		self.name['h5']     = out_dir+'/'+set_['name'].iat[0]+'.h5'

	def is_already_done(self) : 
		''' check if this h5 file is already done/being done'''
		print('')
		dd.dispc('   attempting to work on :'+self.name['h5'],'c','b')
		if os.path.isfile(self.name['h5']) : 
			dd.dispc('     '+self.name['h5']+' already exist','r','r')
			return True 
		if os.path.isdir(self.name['dir']) : 
			dd.dispc('     '+self.name['dir']+' already exist','r','r')
			return True 
		if len(glob.glob(self.name['search']))> 0 :
			dd.dispc('     '+self.name['search']+' already exist','r','r')
			return True 
		return False 


	def create_output_dir_and_lock_file(self) :
		pyos.mkdir(self.name['dir'])
		create_lock_file(self.name['lock'])


	def init(self,inu,fe) :
		''' initialize parameters, matrix that will contain the CC. Do not distinguish btw what depends on this set or not''' 
		# load list of station pair to be correlated from db_cc.h5 
		db_name = self.name['out_dir'] +'/db_cc.h5'
		dset    = 'cc/'+self.name['dir'].split('/')[-2]+'/'+self.name['dir'].split('/')[-1]
		self.cpl = pd.read_hdf(db_name,dset)

		# convenient variable : 
		npath = len(self.cpl) 
		ncmp  = len(inu['cc_cmp'])
		nhr   = len(inu['hr1'])
		nwf   = int(inu['cc_maxlag']*fe)*2+1
						
		# all variables that control the execution of the correlation are in self.ex : 
		self.ex          = {} 
		self.ex['nwf']   = nwf 
		self.ex['npath'] = npath 
		self.ex['ncmp']  = ncmp 
		self.ex['nhr']   = nhr

		# determine list of days that will be correlated it is the same for all path/set:
		self.ex['day_list']=[] 
		for idate_group in inu['date']:
			for idate in idate_group :
				self.ex['day_list'].append((UTCDateTime(idate)).timestamp)

		nday = len(self.ex['day_list'])
		self.ex['ndate_computed'] = nday * nhr 
		self.ex['nday'] = nday 

		# how many date will be saved into the h5 file 
		if inu['save'] == 'day' : 
			self.ex['ndate_save'] = nday 
		if inu['save'] == 'hour' : 
			self.ex['ndate_save'] = self.ex['ndate_computed']
		if inu['save']=='ref' : 
			self.ex['ndate_save'] = 0 

		# Get the list of day at which we will write the result to disk :
		if self.ex['ndate_save'] > 0 :  # means that we will write either daily or hourly cc 
			write_list        = [*range(0,self.ex['nday'],inu['write_step'])]
			self.ex['I_day']  = np.zeros((2,len(write_list)),dtype=np.int)
			self.ex['I_day'][0,:]= np.array(write_list)
			#self.ex['I_day']  = np.array((write_list,write_list))
			self.ex['I_day'][1,:]  = self.ex['I_day'][0,:] + inu['write_step']
			self.ex['I_day'][1,-1] = self.ex['nday']
		else :  # we only save the ref 
			self.ex['I_day']=np.zeros((2,1),dtype='int') 
			self.ex['I_day'][0,0]=0 
			self.ex['I_day'][1,0]=self.ex['nday']

		# initialize ref, ref_nsatck, cc_nstack matrix that will be updated will correlating : 
		self.h5={}
		self.h5['ref'] = np.zeros((ncmp, npath, nwf),dtype=inu['cc_dtype'])
		self.h5['ref_nstack'] = np.zeros((ncmp,npath)) # number of date stacked 
		if inu['save'] != 'ref' :
			self.h5['cc_nstack'] = np.zeros((ncmp,npath,self.ex['ndate_save']))


	def init_h5_file_and_get_mdc_I1(self,db_sta,inu) : 
		# build cc metadata 
		I_sta1 = self.cpl['I_sta1'].to_list()
		I_sta2 = self.cpl['I_sta2'].to_list()
		md={}
		for fd in ['lon','lat','depth','elev','id']  :
			fd1 = db_sta.sta[self.I_set1][fd].iloc[I_sta1].to_numpy()
			fd2 = db_sta.sta[self.I_set2][fd].iloc[I_sta2].to_numpy()
			md[fd]=array((fd1,fd2)).transpose()
		md['id'] = array(md['id'],dtype='S')
		# read md_c : 
		md_c = h5.filename_read_group_as_dict(self.name['out_dir']+'/db_cc.h5','md_c') 

		self.I1 = md_c['I1']
		self.I2 = md_c['I2']

		# write everything to disk :
		dd.dispc('    initializing '+self.name['h5_tmp'],'c','n')
		h5.filename_write_dict_to_group(self.name['h5_tmp'],{'md' :md, 'in_' : inu, 'md_c' : md_c},'/')

		# initializing tree of daily CC
		if inu['save']== 'day' or inu['save']=='hour' :
			ff=h5py.File(self.name['h5_tmp'],'a')
			array_size=(0,self.ex['nwf'])
			array_max_size=(self.ex['ndate_save'],self.ex['nwf'])
			for ipath in range(0,self.ex['npath']) :
				for icmp in inu['cc_cmp'] :
					dset_name ='/cc/'+str(md['id'][ipath,0])[2:-1]+'/'+str(md['id'][ipath,1])[2:-1]+'/'+icmp
					ff.create_dataset(dset_name,array_size,maxshape=array_max_size,dtype=inu['cc_dtype'], chunks=True)
			ff.close()
		

	def open_noise_file(self,kdate,db_sta) :
		''' open noise data (.h5) for this day/network and return h5a, h5b= h5py object'''
		h5a = False 
		h5b = False 
		cdate = UTCDateTime(self.ex['day_list'][kdate]) 
		cyear = str(cdate.year)
		cday  = str(cdate.julday).zfill(3)

		h5_name0 = db_sta.set['path'][self.I_set1]+'/'+cyear+'/day_'+cday+'.h5'
		h5_name1 = db_sta.set['path'][self.I_set2]+'/'+cyear+'/day_'+cday+'.h5'
		if not os.path.isfile(h5_name0) : 
			dd.dispc('     '+h5_name0+'  not found','r','n')
			return h5a, h5b
		if  not os.path.isfile(h5_name1) : 
			dd.dispc('      '+h5_name1+'  not found','r','n')
			return h5a, h5b 
		try :
			h5a = h5py.File(h5_name0,'r') 
			h5b = h5py.File(h5_name1,'r') 
		except :
			return h5a, h5b 
		return h5a, h5b 


	def read_data(self,h5a,h5b,kcmp,kcpl,db_sta) : 
		''' read only noise data that will be correlated from an open h5_file  '''

		trace0 = []
		trace1 = [] 

		# get indice of the set and stations : 
		I_set0 = self.I_set1 
		I_set1 = self.I_set2 
		I_sta0 = self.cpl.loc[kcpl,'I_sta1']
		I_sta1 = self.cpl.loc[kcpl,'I_sta2']

		# now we can get the station name from db_sta :
		sta_name0 = db_sta.sta[I_set0].loc[I_sta0,'id']
		sta_name1 = db_sta.sta[I_set1].loc[I_sta1,'id']

		# get the dset name : 
		net0 = sta_name0.split('.')[0]
		net1 = sta_name1.split('.')[0]
		dset0 = net0+'/'+sta_name0.split('.')[1]+'.'+sta_name0.split('.') [2]+'/'+kcmp[0]
		dset1 = net1+'/'+sta_name1.split('.')[1]+'.'+sta_name1.split('.') [2]+'/'+kcmp[1]

		# attempt to read the trace  :
		try :
			if not dset0 in h5a : 
				return trace0, trace1 	
			if not dset1 in h5b  :
				return trace0, trace1
		except : 
			return trace0, trace1

		try :
			trace0 = h5a[dset0][:]
			trace1 = h5b[dset1][:]
		except :
			pass 
		return trace0, trace1 

	def write_cc(self,db_sta,cmp_list,I1,I2) :
		''' at the end of a group of day : write the daily/hourly cc'''
		ff = h5py.File(self.name['h5_tmp'],'a')
		for ipath in range(0,self.ex['npath']) :
			kcmp=-1
			for icmp in cmp_list:
				kcmp=kcmp+1
				I_sta1 = self.cpl.loc[ipath,'I_sta1']
				I_sta2 = self.cpl.loc[ipath,'I_sta2']
				id1 = db_sta.sta[self.I_set1].loc[I_sta1,'id']
				id2 = db_sta.sta[self.I_set2].loc[I_sta2,'id']
				dset_name = '/cc/'+id1+'/'+id2
				dset_name = dset_name + '/'+icmp
				#dset_name ='/cc/'+str(self.cc_info['id'][ipath,0])[2:-1]+'/'+str(self.cc_info['id'][ipath,1])[2:-1]
				#dset_name = dset_name + '/'+icmp
				ff[dset_name].resize((I2,self.ex['nwf']))
				ff[dset_name][I1:I2,:] = self.h5['cc'][kcmp,ipath,:,:]
		ff.close()


	def write_all(self,db_sta, cmp_list) : 
		''' once all days are correlated : write the reference, cc_nstack, ref_nstack'''
		ff = h5py.File(self.name['h5_tmp'],'a')
		for ipath in range(0,self.ex['npath']) :
			kcmp =-1
			for icmp in cmp_list:
				kcmp = kcmp+1 
				I_sta1 = self.cpl.loc[ipath,'I_sta1']
				I_sta2 = self.cpl.loc[ipath,'I_sta2']
				id1 = db_sta.sta[self.I_set1].loc[I_sta1,'id']
				id2 = db_sta.sta[self.I_set2].loc[I_sta2,'id']
				dset_name = id1+'/'+id2 
				#dset_name =str(self.cc_info['id'][ipath,0])[2:-1]+'/'+str(self.cc_info['id'][ipath,1])[2:-1]
				dset_name = dset_name + '/'+icmp
				ff['/ref/'+dset_name]=self.h5['ref'][kcmp,ipath,:]/ max(1,self.h5['ref_nstack'][kcmp,ipath])
		ff['/ref_nstack']=self.h5['ref_nstack'] 
		if 'cc_nstack' in self.h5 :
			ff['/cc_nstack'] = self.h5['cc_nstack']
		ff.close()

	def clean_all(self) : 
		''' at the end : remove lock file, temoporary dir, and move the h5 file to its final name'''
		pyos.rename(self.name['h5_tmp'],self.name['h5'])
		pyos.remove(self.name['lock'])
		pyos.remove_dir(self.name['dir'])



#============================================================================
class DB_sta : 
	def __init__(self,out_dir) :
		''' import directly fe, cut_len, set, year, and add manually sta which is a list of dataframe'''
		db_name = out_dir+'/db_sta.h5'
		h5.import_h5_data_to_class(self,db_name)
		self.sta=[]
		for iset in self.set['name'] :
			self.sta.append(pd.read_hdf(db_name,'sta/'+iset))
		#self.sta = pd.concat(self.sta,ignore_index=True)




##########################################################################################
##########################################################################################
###########         ####            ####           ######           ######################
##########          ####            ####            #####            #####################
##########     #########    ####    ####    ####    #####    ####    #####################
#########     ##########    ####    ####    ####   ######    ####   ######################
#########     ##########    ####    ####           ######           ######################
#########      #########    ####    ####    ####    #####    ####    #####################
##########          ####            ####    ####    #####    ####    #####################
###########         ####            ####    #####    ####    #####    ####################
###################################/#######################################################
##########################################################################################


#
def correlate_this_path(trace1,trace2,maxlag,I1,I2,cc_type,data_type) : 
    ''' correlate a single path at a give date. Low lvl functions only'''
    fh = globals()[cc_type]
    trace1[np.isnan(trace1)] = 0.0
    trace2[np.isnan(trace2)] = 0.0
    trace01s = [trace1[I1[itr]:I2[itr]]  for itr in range(0,len(I1),1)]
    trace02s = [trace2[I1[itr]:I2[itr]]  for itr in range(0,len(I1),1)]
    maxlag   = int(maxlag)
    N        = len(trace01s)
    corr     = np.zeros([N,2*maxlag+1], dtype=data_type)
    ncorr    = np.zeros(N)
    for ihr in range(N):
        if trace01s[ihr].std() > 1e-15:      #  numpy.max(numpy.abs(trace01s[i])) >= 1e-20:
            if trace02s[ihr].std() > 1e-15:  #  numpy.max(numpy.abs(trace02s[i])) >= 1e-20:
                #try:
                corr[ihr] = fh(trace02s[ihr], trace01s[ihr],maxlag).astype(data_type)
                ncorr[ihr] = 1 # TO BE MODIFIED !!! ncorr[ihr] = 1 or 0
                #except:pass
    return corr,ncorr



 
def ctp_xcorr_norm(trace01, trace02,maxlag): # cross-coherence Time normalized before Corr
    lentrace   = len(trace01)
    goodnumber = int(maxlag+lentrace)
    tr2 = np.zeros(goodnumber)
    tr2[0:lentrace] = trace02
    tr2[0:lentrace] /= np.sqrt(np.sum(tr2[0:lentrace]**2))
    tr2 = fftpack.fft(tr2,overwrite_x=True)
    tr2.imag *= -1
    tr1 = np.zeros(goodnumber)
    tr1[maxlag:int(maxlag+lentrace)]= trace01
    tr1[maxlag:int(maxlag+lentrace)] /= np.sqrt(np.sum(tr1[maxlag:int(maxlag+lentrace)]**2))
    tr2 *= fftpack.fft(tr1,overwrite_x=True)
    tr2[np.isnan(tr2)] = 0.0+0.0j;
    tr2[np.isinf(tr2)] = 0.0+0.0j;
    return (fftpack.ifft(tr2,overwrite_x=True)[0:2*maxlag+1].real)

def ctp_xcorr(trace01, trace02,maxlag): # cross-correlation
    lentrace=len(trace01)
    goodnumber = int(maxlag+lentrace)
    tr2 = np.zeros(goodnumber)
    tr2[0:lentrace] = trace02
    tr2 = fftpack.fft(tr2,overwrite_x=True)
    tr2.imag *= -1
    tr1 = np.zeros(goodnumber)
    tr1[maxlag:int(maxlag+lentrace)]= trace01
    tr2 *= fftpack.fft(tr1,overwrite_x=True)
    tr2[np.isnan(tr2)] = 0.0+0.0j;
    tr2[np.isinf(tr2)] = 0.0+0.0j;
    return (fftpack.ifft(tr2,overwrite_x=True)[0:2*maxlag+1].real)

def ctp_coher(trace01, trace02,maxlag): # cross-coherence
    lentrace   = len(trace01)
    goodnumber = int(maxlag+lentrace)
    tr2 = np.zeros(goodnumber)
    tr2[0:lentrace] = trace02
    tr2 = fftpack.fft(tr2,overwrite_x=True)
    tr2.imag *= -1
    tr1 = np.zeros(goodnumber)
    tr1[maxlag:int(maxlag+lentrace)]= trace01
    tr1 = fftpack.fft(tr1,overwrite_x=True)
    ## wtlvl
        #wlvl = 1.0
        #temp2 = numpy.absolute(tr2)
        #temp2[numpy.isnan(temp2)] = 0.0;
        #temp2[numpy.isinf(temp2)] = 0.0;
        #mtemp2 = numpy.amax(temp2) * wlvl / 100
        #numpy.putmask(temp2,temp2 <= mtemp2, mtemp2)
        #temp1 = numpy.absolute(tr1)
        #temp1[numpy.isnan(temp1)] = 0.0;
        #temp1[numpy.isinf(temp1)] = 0.0;
        #mtemp1 = numpy.amax(temp1) * wlvl / 100
        #numpy.putmask(temp1,temp1 <= mtemp1, mtemp1)
    ##
    tr2 = (tr1 * tr2) / (np.absolute(tr1) * np.absolute(tr2))
    #tr2 = (tr1 * tr2) / (temp1 * temp2)
    tr2[np.isnan(tr2)] = 0.0+0.0j;
    tr2[np.isinf(tr2)] = 0.0+0.0j;
    return (fftpack.ifft(tr2,overwrite_x=True)[0:2*maxlag+1].real)

def ctp_coher_tap(trace01, trace02,maxlag): # cross-coherence + taper
    taper      = 0.15
    strength   = 1000
    lentrace   = len(trace01)
    goodnumber = int(maxlag+lentrace)
    tr2 = np.zeros(goodnumber)
    tr2[0:lentrace] = trace02
    tr2 = fftpack.fft(tr2,overwrite_x=True)
    tr2.imag *= -1
    tr1 = np.zeros(goodnumber)
    tr1[maxlag:maxlag+lentrace]= trace01
    tr1 = fftpack.fft(tr1,overwrite_x=True)
    tt = np.append(signal.tukey(divmod(len(tr1),2)[0],taper),signal.tukey(divmod(len(tr1),2)[0] + divmod(len(tr1),2)[1],taper))
    tt = 1+(1-tt)*strength
    tr2 = (tr1 * tr2) / (np.absolute(tr1) * np.absolute(tr2) *tt)
    tr2[np.isnan(tr2)] = 0.0+0.0j;
    tr2[np.isinf(tr2)] = 0.0+0.0j;
    return (fftpack.ifft(tr2,overwrite_x=True)[0:2*maxlag+1].real)

def ctp_deconv(trace01, trace02, maxlag): # deconvolution
    ''' voir la signification d'EXTRA!!!'''
    EXTRA      = 1
    lentrace   = len(trace01)
    goodnumber = int(maxlag+lentrace)
    tr2 = np.zeros(goodnumber)
    tr2[0:lentrace] = trace02
    tr2 = fftpack.fft(tr2,overwrite_x=True)
    tr2.imag *= -1
    tr1 = np.zeros(goodnumber)
    tr1[maxlag:maxlag+lentrace]= trace01
    tr1 = fftpack.fft(tr1,overwrite_x=True)
    ################################################################
    ##### WATER LVL
    #
    ############
    #wlvl = EXTRA
    ############
    #temp = numpy.absolute(tr2)
    #mtemp = numpy.amax(temp)* wlvl / 100
    #numpy.putmask(temp,temp <= mtemp, mtemp)
    #tr2 = (tr1 * tr2) / (temp ** 2)
    ################################################################
    ##### NOISE PADDING
    #tr2_n = (numpy.random.random(goodnumber)-0.5)*2*(numpy.std(trace02)*EXTRA/100)
    tr2_n = (np.random.random(goodnumber)-0.5)*2*(EXTRA)
    tr2_n[0:lentrace] = trace02
    tr2_n = fftpack.fft(tr2_n,overwrite_x=True)
    tr2_n.imag *= -1
    tr2 = (tr1 * tr2) / (np.absolute(tr2_n)**2)
    #################################################################
    #################################################################
    ##### SMOOTHING
    #temp = numpy.absolute(tr2)
    ###########
    #window_len = EXTRA
    ###########
    #s = numpy.r_[2*temp[0]-temp[window_len-1::-1],temp,2*temp[-1]-temp[-1:-window_len:-1]]
    ##s = numpy.r_[numpy.zeros(len(temp[window_len-1::-1])),temp,numpy.zeros(len(temp[-1:-window_len:-1]))]
    #w = numpy.ones(window_len,'d')
    #y = numpy.convolve(w/numpy.sum(w),s,mode='same')
    #temp = y[window_len:-window_len+1]
    #tr2 = (tr1 * tr2) / (temp**2)
    #################################################################
    #################################################################
    return (fftpack.ifft(tr2,overwrite_x=True)[0:2*maxlag+1].real)


def ctp_xcorr_norm_envelope(trace01, trace02,maxlag): 
    # cross-coherence Time normalized before Corr
    # since the 'envelope' is always larger than zero we need add a demean to make the correlation meaningful
    # THIS IS TOO SIMPLE TO BE EXPLAINED!!!!
    trace01 = signal.detrend(trace01, type='constant')
    trace02 = signal.detrend(trace02, type='constant')
    lentrace   = len(trace01)
    goodnumber = int(maxlag+lentrace)
    tr2 = np.zeros(goodnumber)
    tr2[0:lentrace] = trace02
    tr2[0:lentrace] /= np.sqrt(np.sum(tr2[0:lentrace]**2))
    tr2 = fftpack.fft(tr2,overwrite_x=True)
    tr2.imag *= -1
    tr1 = np.zeros(goodnumber)
    tr1[maxlag:int(maxlag+lentrace)]= trace01
    tr1[maxlag:int(maxlag+lentrace)] /= np.sqrt(np.sum(tr1[maxlag:int(maxlag+lentrace)]**2))
    tr2 *= fftpack.fft(tr1,overwrite_x=True)
    c12 = (fftpack.ifft(tr2,overwrite_x=True)[0:2*maxlag+1].real)
    return c12






#=========================================================================

def get_file_indice(out_dir) : 
	# get the number of files to be correlated from db_cc.h5
	ff = h5py.File(out_dir+'/db_cc.h5','r')
	nfile = ff['nset'][()]
	ff.close()
	# get a randomize file_list indice : 
	I_file = [*range(0,nfile)]
	random.shuffle(I_file)
	return I_file

def get_and_create_directory_name(in_) :
	''' determine the name of the output directory C1_* from user input'''
	#get script that called the code name :
	script_name = inspect.stack()[2]
	caller = script_name[1][0:-3].split('/')[-1]
	if 'out_dir' in in_ :
		out_dir = in_['out_dir']+'/'
	else :
		out_dir ='./'
	out_dir  = out_dir+in_['path_out']+'C1_'+caller+'_'+ in_['tag'] 
	out_dir +='__'+in_['cc_func'][4:-1]+in_['cc_func'][-1]  

	# create the C1 dir : 
	if not os.path.isdir(out_dir) : 
			pyos.mkdir(out_dir)

	dd.dispc('','y','b')
	dd.dispc('-----------------------------------------------------','y','b')
	dd.dispc('  correlations will be saved in '+str(out_dir),'y','n')
	return out_dir 

def create_lock_file(filename) : 
    ff=open(filename,'wb')
    pickle.dump(filename,ff)
    ff.close()


	
	
