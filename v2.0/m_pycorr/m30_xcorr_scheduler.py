################################################
# m30_xcorr_scheduler.py 
# Laurent Stehly (UGA)
# Pierre Boue (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################

import time,os, inspect, sys, h5py, random, glob, datetime, pickle, ipdb 
import pandas as pd, numpy as np
from obspy.core import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from m_pycorr.mods import h5, lang, dd, pyos 
import ipdb

# Ne supporte les sets commencant par des chiffre => les renommer en lettre ?

#-------------------------------------------------------------------------------------------------------------------
# The role of the scheduler is to create the following files :
# - C1/db_sta.h5 : panda tables containing station metadata built from data_1.0hz/daily/net/year/db.pkl files 
# - C1/db_cc.h5  : list of indice of set and stations to be correlated file per file + all other info to computed CC
# - C1/db_cc_file : used to read correlations : indice of station pair in each file (not implemented now)
#
#-------------------------------------------------------------------------------------------------------------------
#             LOGIC  FLOW 
#---------------------------------------------
# if db_sta.h5 does not exist : we create it 
#              do exist : we do nothing, and do not attempt to read it
#              is locked : wait for it to exist
#
# Once db_sta.h5 exist (has been created by one of the process) we continue : 
#
# if db_cc.h5 : does not exist : we read db_sta.h5 if we did not create it ourself, and build db_cc.h5 
#             : does exist : we do nothing, and do not attempt to read it 
#             : is locked : we wait for it to exist 
#
# Once we are sure that both db_sta.h5 and db_cc.h5 exist, we generate others db_cc_... file and go to xcorr 
# db_sta.init and db_cc.init are set to True or False depending if db_sta.h5/db_cc.h5 has been read/constructed or are # empty. So we use it to know we should read db_sta.h5/db_cc.h5 or not 
#--------------------------------------------------------------------------------------------------------------------


# TODO : modidier les tables pour avoir l'id/nom de fichier comme index ?

def schedule_all(in_,out_dir) :

	db_sta = DB_sta(out_dir,in_['path'],in_['date']) 
	if db_sta.init : 
		db_sta.disp_station_info() 
		db_sta.disp_station_list(in_['verbose_lvl']) 

	while not os.path.isfile(db_sta.name['h5']) : 
		pass 
#	db_sta.read()#
#	db_sta.sta[0]
	db_cc = DB_cc(out_dir,db_sta,in_) 
	if db_cc.init : 
		db_cc.disp_path_repartition_per_file(in_['dist'],in_['verbose_lvl'])

	while not os.path.isfile(db_cc.name['h5']) :
		pass 

#	generate_db_c1(db_sta,db_cc,in_,out_dir)
#	generate_db_c1_file(db_sta,db_cc,out_dir) 

	if in_['verbose_lvl'] >= 2 :
		print('\n \n')
		dd.main(db_cc,lvl=in_['verbose_lvl']-1)
		print('\n \n')
		dd.main(db_sta,lvl=in_['verbose_lvl']-1)

	if in_['verbose_lvl'] ==4 :
		db_sta.read()
		db_sta.disp_full_station_list()


#################################################################################################################
#
#
#  GENERATE C1.../db_file.h5 files 
#
#
#############################################################################################
def generate_db_c1(db_sta,db_cc,in_,out_dir) :
	db_name = out_dir +'/db_c1.h5' 
	db  = _init_name(out_dir,'db_c1')

	# does this db already exist ?
	if len(glob.glob(db['name']['search'])) > 0: 
		return 

	# no so lets go 
	create_lock_file(db['name']['lock'])

	if db_sta.init == False : 
		db_sta.read()
	if db_cc.init == False :
		db_cc.read()

	# create a table whose index = filename and a single columns which is the number of path stored in this file
	nset = len(db_cc.cc)
	ncpl = np.zeros((nset),dtype='int32')
	for iset in range(0,nset) :	
		ncpl[iset] = len(db_cc.cc[iset]) 

	xcorr_file_table = pd.DataFrame(ncpl,index=db_cc.set['name'],columns=['ncpl'])
	xcorr_file_table.to_hdf(db['name']['h5'],'file',format='table')

	# create a table whose index =id of the station : 
	sta_table = pd.concat(db_sta.sta,ignore_index=True)
	sta_table.index = sta_table['id']
	sta_table = sta_table.reindex(['lon','lat','elev','depth'],axis='columns')
	sta_table.to_hdf(db['name']['h5'],'sta',format='table')
	
	# list of cmp 
	pd.DataFrame(in_['cc_cmp']).to_hdf(db['name']['h5'],'cmp',format='table')
	pyos.remove(db['name']['lock'])


def generate_db_c1_file(db_sta,db_cc,out_dir) : 
	db_name = out_dir +'/db_c1_by_file.h5' 
	db = _init_name(out_dir,'db_c1_by_file')

	# does this db already exist ?
	if len(glob.glob(db['name']['search'])) > 0 : 
		return 

	#no so lets go : 
	create_lock_file(db['name']['lock']) 
	if db_sta.init == False :
		db_sta.read()
	if db_cc.init==False :
		db_cc.read()
	
	J_sta = db_sta._get_indice_first_station_in_each_set() 

	# loop on files i.e db_cc_file.h5 group :
	ff = h5py.File(db['name']['h5'],'w')
	h5_sta   =[] # 2D List of station indice that will be stored into the h5 file :
	h5_dist = [] # 1D list of distance 
	h5_az   = [] 
	h5_baz  = [] 
	# loop on path that will be stored in each file :
	for ifile, vfile in db_cc.set.iterrows() :
		group_name = vfile['name']
		I_set1 = vfile['I_set1'] 
		I_set2 = vfile['I_set2']  
		J_sta1 = J_sta[I_set1]   # will be used to convert relative to absolute indice of sta1
		J_sta2 = J_sta[I_set2]   # will be used  to convert relative to absolute indice of sta2
		# loop on station couple, need to get their absolute indice ... 
	
		h5_sta = np.array((db_cc.cc[ifile]['I_sta1'].array+J_sta1,db_cc.cc[ifile]['I_sta2'].array+J_sta2)).transpose()

		for icpl in range(0,len(h5_sta)) :
			dist_,az,baz = db_sta.get_dist_and_azimuth(I_set1,h5_sta[icpl,0]-J_sta1,I_set2,h5_sta[icpl,1]-J_sta2)	
			# convertir indice relatif en absolu + stockage en list 1D/2D
			h5_dist.append(dist_)
			h5_baz.append(baz)
			h5_az.append(az)
		# write them on disk :

		ff.create_dataset(group_name+'/I_sta',data=h5_sta)
		ff.create_dataset(group_name+'/az',data=h5_az)
		ff.create_dataset(group_name+'/baz',data=h5_baz)
		ff.create_dataset(group_name+'/dist',data=h5_dist)
	ff.close()
	pyos.remove(db['name']['lock'])


def _init_name(out_dir,db_name) :
	db={}
	db['name']={}
	db['name']['h5']     = out_dir+'/'+db_name+'.h5' 
	db['name']['lock']   = out_dir+'/'+db_name+'.'+str(os.getpid())+'.lock' 
	db['name']['search'] = out_dir+'/'+db_name+'.*'
	return db  



#######################################################################################################################
#
#  DB CC
#   
########################################################################################################################
class DB_cc :
	def __init__(self,out_dir,db_sta,in_): 
		self.name = {}
		self.name['h5'] = out_dir+'/db_cc.h5'
		self.name['lock'] = out_dir+'/db_cc.'+str(os.getpid())+'.lock'
		self.name['search'] = out_dir+'/db_cc.*'
		self.init = False 

		# db_cc.h5 file already exist => return 
		if os.path.isfile(self.name['h5']): 
			return 

		#if lock file does not exist => create it : means we will write to disk 
		if len(glob.glob(self.name['search'])) == 0 :
			create_lock_file(self.name['lock'])
			dd.dispc('  '+self.name['search']+' does not exist => create it','c','n')
			if db_sta.init == False :  
				db_sta.read()
			self.determine_npath_per_h5_file(in_,db_sta.fe) 
			self.get_metadata_common_to_all_cc(in_,db_sta.fe,db_sta.cut_len) 
			self.organize_correlations(db_sta,in_)
			self.init = True 
			self.save()
		

	def save(self) : 
		''' save the data contained into self if did create the lock file '''
		# if the db_cc file already exist we return 
		if os.path.isfile(self.name['h5']) :
			dd.dispc('  db_cc.h5 already exist => returning','r','r')
			return 
		# did we create the lock file ? if so save. Otherwise do nothing 
		if os.path.isfile(self.name['lock']):
			dd.dispc('  saving db_cc.h5','r','r')
			h5.save_class_data_to_h5(self,self.name['h5']) 
			# need to save manually list of dataframe 
			for icc in range(0,len(self.cc)) : 
				self.cc[icc].to_hdf(self.name['h5'],'cc/'+self.set['name'][icc])
			pyos.remove(self.name['lock'])


	def read(self) : 
		h5.import_h5_data_to_class(self,self.name['h5']) 
		# read manually cc which is a list of dataframe. Group and subgroup form the name of
		self.cc=[]
		for ifile in self.set['name'] : 
			self.cc.append(pd.read_hdf(self.name['h5'],'/cc/'+ifile))


	def determine_npath_per_h5_file(self,in_,fe) :
		''' determine how many path will be stored in a .h5 file depending on user input.'''

		# size of a single path, all cmp, single date in bytes (not Go)
		size_single_path_date = (in_['cc_maxlag']*2.0+1.0) * fe * len(in_['cc_cmp'])
		if in_['cc_dtype']=='float64'    : size_single_path_date = size_single_path_date * 8.0
		elif in_['cc_dtype']=='float32'  : size_single_path_date = size_single_path_date * 4.0
		elif in_['cc_dtype']=='float16'  : size_single_path_date = size_single_path_date * 2.0
		else : 
			dd.dispc(' cc_dtype different from float16/32/64','y','r')
		
		# get number of hour per day and the number of days
		nhr  = len(in_['hr1'])
		nday = 0
		for idate_group in in_['date'] :
			nday = nday + len(idate_group)

		if in_['save'] == 'ref' :
			ndate = 1 
		elif  in_['save'] == 'day' :
			ndate = nday + 1
		elif  in_['save'] == 'hour' :
			ndate = nday*nhr + 1 

		# size in Go of a single path for all date and all component :
		size_single_path=size_single_path_date*ndate/1024./1024./1024. 
		npath_per_file=int(in_['file_size']/size_single_path)
		
		# we store all this kind of information in ex :
		self.ex= {}
		self.ex['npath_per_file'] = max(1,npath_per_file )
		self.ex['size_single_path'] = size_single_path 

		dd.dispc('','y','b')
		dd.dispc('-------------- npath per file -----------------------','y','b')
		dd.dispc('  we have '+str(ndate)+' dates to save and '+str(len(in_['cc_cmp']))+' cmp','y','n')
		dd.dispc('  => You asked to have files of at most '+str(in_['file_size'])+ 'Go','y','n')
		dd.dispc('  => we will put at most '+str(self.ex['npath_per_file'])+' paths per h5 file','y','n')
		dd.dispc('    = '+str(round(self.ex['npath_per_file']*size_single_path,2))+' Go per file','y','n')
		dd.dispc('  note that it is possible to have a file size < to the size of a single path','y','n')
		dd.dispc('  in that case we will exceed the filesize that you asked','y','n')

		if in_['save']=='ref' :
			size_ram = size_single_path_date * self.ex['npath_per_file'] /1024/1204/1024. 
		if in_['save']=='day' :
			size_ram = size_single_path_date * (in_['write_step']+1)*self.ex['npath_per_file'] /1024/1204/1024. 
		if in_['save']=='hour' :
			size_ram = size_single_path_date * (in_['write_step']*nhr+1)*self.ex['npath_per_file'] /1024/1204/1024. 

		self.ex['size_ram'] = size_ram 
		dd.dispc('  MEMORY USAGE : about '+ str(round(size_ram,2))+' Go','y','b')


	def get_metadata_common_to_all_cc(self,in_,fe,cut_len) :
		'''set self.md_c = metadata common to all CCs determined from user input, fe and cut_len '''
		md_c={}
		md_c['tau']  =1./fe
		# information on noise data :
		time={}
		time['fe']   = fe 
		time['t1']   = cut_len
		time['t2']   = 86400-time['t1']-1             
		time['npts'] = time['fe']*(86400-time['t1']-time['t1'])        #
		time['time'] = np.arange(time['t1'],time['t2'],1.0/time['fe']) # in second w.r to midnight
		ds1=np.array(in_['hr1'])*3600-time['t1']                       # difference de temps [s]
		ds2=np.array(in_['hr2'])*3600-time['t1']                       # entre chaque fen et le t1 des traces

		#keep the indice of each hourly slice :
		md_c['I1']=np.round(ds1*time['fe']).astype('int')# idem ms en nb de points
		md_c['I2']=np.round(ds2*time['fe']).astype('int')
		if md_c['I2'][-1] > time['npts'] :
			dd.dispc('  sch_get_metadata : will try to correlate time that do not exist','r','b')

		md_c['version']= 0.9
		md_c['t'] = np.linspace(-in_['cc_maxlag'],in_['cc_maxlag'],int(2*in_['cc_maxlag']*1/md_c['tau']+1))

		#add the list of components from in_ : so that it can be changed later if we rotate the CC : 
		md_c['cmp']=in_['cc_cmp']

		#construct date vector = list of day correlated : 
		md_c['date1']=[]
		md_c['date2']=[]

		if in_['save'] == 'day' or in_['save']=='ref':
			for idate_group in in_['date']:
				for idate in idate_group :
					md_c['date1'].append((UTCDateTime(idate)+time['t1']).timestamp)
					md_c['date2'].append((UTCDateTime(idate)+time['t2']).timestamp)
	
		if in_['save'] =='hour' :
			for idate_group in in_['date'] :
				for idate in idate_group :
					for ihr in range(0,len(in_['hr1'])) :
						md_c['date1'].append((UTCDateTime(idate)+(in_['hr1'][ihr]*3600)).timestamp)
						md_c['date2'].append((UTCDateTime(idate)+(in_['hr2'][ihr]*3600)).timestamp)
	
		#add a matlab - readable version :
		md_c['date1_ordinal']=[]
		md_c['date2_ordinal']=[]
		for idate in range(0, len(md_c['date1'])) : 
			cdate=UTCDateTime(md_c['date1'][idate]) 
			mdate1=cdate.toordinal()+(cdate.hour*3600+cdate.minute*60+cdate.second)/86400.
			cdate=UTCDateTime(md_c['date2'][idate]) 
			mdate2=cdate.toordinal()+(cdate.hour*3600+cdate.minute*60+cdate.second)/86400.
			md_c['date1_ordinal'].append(mdate1)
			md_c['date2_ordinal'].append(mdate2)
		self.md_c = md_c 


	def organize_correlations(self,db_sta,inu) :
		''' set self.cc and self.set : determine which couple of station will be stored in which file'''
		
		self.cc = [] # 1D list of data frame containing the index of the station 
		set_list =[]
		npath_per_file = self.ex['npath_per_file']

		nset = len(db_sta.set) 
		npath_wrong_dist = 0 
		# correlations croisees : (FR-CH, FR-IV,....)
		if inu['outer'] :
			for iset1 in range(0,nset) :
				for iset2 in range(iset1+1,nset) :
					if not db_sta.is_there_a_common_year(iset1,iset2) : 
						continue 
					name1 = db_sta.set.loc[iset1,'name']
					name2 = db_sta.set.loc[iset2,'name']
#					print('  db_cc : constructing '+name1+'  '+name2)
					db_cc1, nwrong_dist = self._construct_db_cc(db_sta,iset1,iset2,inu,npath_per_file,'outer') 
					if len(db_cc1) > 0 :
						set_list = self._append_cc_list(db_cc1,iset1,iset2,set_list,name1,name2,'xcorr_all/xcorr') 

	    # correlation interne a un set (FR-FR, CH-CH,...)
		if inu['inner'] :
			for iset1 in range(0,nset) :
				name1 = db_sta.set.loc[iset1,'name']
				name2 = name1 #db_sta.set.loc[iset1,'name']
				db_cc1, nwrong_dist = self._construct_db_cc(db_sta,iset1,iset1,inu,npath_per_file,'inner')
				if len(db_cc1) > 0 :
					npath_wrong_dist = npath_wrong_dist + nwrong_dist
					set_list = self._append_cc_list(db_cc1,iset1,iset1,set_list,name1,name2,'xcorr_all/xcorr') 

		# a set vs others set of stations 
		if len(inu['set_vs_all']) > 0 :
			for jset1 in inu['set_vs_all'] :
				iset1 = np.where(db_sta.set['name'].isin([jset1])==True)[0][0]
				for iset2 in range(0,nset) : 
					if iset1 == iset2 :
						continue 
					name1 = db_sta.set.loc[iset1,'name']
					name2 = db_sta.set.loc[iset2,'name']
					db_cc1, nwrong_dist = self._construct_db_cc(db_sta,iset1,iset2,inu,npath_per_file,'outer')
					if len(db_cc1) > 0 :
						npath_wrong_dist = npath_wrong_dist + nwrong_dist
						set_list = self._append_cc_list(db_cc1,iset1,iset2,set_list,name1,name2,'xcorr_net/net') 
	
		# CC between a station vs stations belongin to others set of stations 
		# we first loop on set defined in in_['sch_vs_net']
		kset =-1 
		for jset in inu['vs_set'] :
			kset=kset+1
			iset1 = np.where(db_sta.set['name'].isin([jset])==True)[0][0] # numero du set ou se trouve la vs 
			vs = inu['vs_list'][kset]                                   # indice des source virtuelles
			for ivs in vs : 
				for iset2 in range(0,nset) : 
					if inu['vs_inner']==False and (iset1 == iset2) : 
						continue 
					name1 = '' #db_sta.set.loc[iset1,'name']
					name2 = db_sta.set.loc[iset2,'name']
					db_cc1, wrong_dist = self._construct_db_cc(db_sta,iset1,iset2,inu,npath_per_file,'outer',sta1_list=[ivs])
					if len(db_cc1) > 0 :
						npath_wrong_dist += wrong_dist
						prefix = 'xcorr_vs/vs_'+db_sta.sta[iset1].at[ivs,'id'].replace('.','_')
						set_list = self._append_cc_list(db_cc1,iset1,iset2,set_list,name1,name2,prefix) 


		self.set = pd.DataFrame(set_list,columns=['I_set1','I_set2','name']) 
		self.nset = len(self.set)


	def _append_cc_list(self,db_cc1,kset1,kset2,set_list,set1_name,set2_name,prefix) :
		''' update self.I_sta and set_list which is a bit weird :-) '''
		knum =-1 
		for ifile in range(0,len(db_cc1)) : 
			knum = knum + 1
			self.cc.append(pd.DataFrame(db_cc1[ifile]['cpl'],columns=['I_sta1','I_sta2']))
			set_list.append([kset1,kset2,prefix+'_'+set1_name+'_'+set2_name+'_'+str(knum).zfill(4)])
		return set_list 


	def _construct_db_cc(self,db_sta,set1,set2,inu,npath_per_file,type_,sta1_list=None) :
		''' private function of organize_cc : for a pair of set of stations, list all CCs that will be computed'''
  		#  For a pair of set of stations return db_cc a 2D list of couple of  stations to be correlated. 
  		#    	- each row, is a 1D list will be stored in a different h5 file. 
	    #       - each element is a dictonnary {'set_indice': [0, 1], 'cpl': [[0, 1], [2, 1], [4, 1]]}
	    #        - the list is 2D since we impose that at most npath_per_file couple of station per h5 file. 
		#	does not really depends on the bd_cc class so it is quite outside the class actually :) 
		#
		# if type_ = 'inner' : we list the couple of stations within a network, i.e the correlation 
		#                     to be computed between FR-FR for instance. In that case we list the couple 
		#					  [0-1,0-2,1-1,1-2,2-2] but not [1-0,2-0,2-1]
		# if type_ = 'outer' : we list the couple of stations between two different network for instance FR-CH, 
		#                      in that case we list everything : [0-1, 0-2, 1-0,1-1,1-2,2-0,2-1,2-2]  
		# if sta1_list is set : that for computing the CC between a virtual source and the rest of the network. 
		#                     : in that case we generate a 1D list per station indice listed in sta1_list. 
		#
		# return an empty list if no couple of station is selected. That happens for instance if there is no station
		# separated by the right distance.
		dist = inu['dist']
		db_cc = []
		kwrong_dist = 0 
		ncpl = -1 
		# loop on station pair 
		nsta1 = len(db_sta.sta[set1])
		nsta2 = len(db_sta.sta[set2])
		# sta1_list = None mean we loop on station of the 1st set
		# sta1_list = one index = we use this station as a virtual source vs the rest of the stations
		if sta1_list == None :
			sta1_list = list(range(0,nsta1))
		for ista1 in sta1_list : # range(0,nsta1) :
			if type_ == 'inner' :
				sta2 = range(ista1,nsta2)
			elif type_ == 'outer' :
				sta2 = range(0,nsta2)	
			for ista2 in sta2 :
				# check distance range + coordinante, otherwise skip this path :
				if dist != None :
					#dist_,az,baz,lon1,lat1,lon2,lat2 = db_sta.get_dist_and_azimuth(set1,ista1,set2,ista2)[0]
					sta_info=db_sta.get_dist_and_azimuth(set1,ista1,set2,ista2)
					if (sta_info[3] < inu['lon'][0]) or (sta_info[5] < inu['lon'][0]) :
						kwrong_dist = kwrong_dist+1 
						continue 
					if (sta_info[3] > inu['lon'][1]) or (sta_info[5] > inu['lon'][1]) :
						kwrong_dist = kwrong_dist+1 
						continue 
					if (sta_info[4] < inu['lat'][0]) or (sta_info[6] < inu['lat'][0]) :
						kwrong_dist = kwrong_dist+1 
						continue 
					if (sta_info[4] > inu['lat'][1]) or (sta_info[6] > inu['lat'][1]) :
						kwrong_dist = kwrong_dist+1 
						continue 
					dist_ = sta_info[0]
					if dist_ < dist[0] or dist_ > dist[1] :
						kwrong_dist= kwrong_dist + 1 
						continue
					# we start  => we init db_cc 
					if ncpl == -1 :
						db_cc= [] 	
						db_cc.append({})
						db_cc[-1]['cpl']=[]
						ncpl = 0 
					# check that we do not put more than npath_per_file couple of station per h5 file
					if ncpl >= npath_per_file  :
						db_cc.append({})
						db_cc[-1]['cpl']=[]
						ncpl=0
					# Add this path to the list of path to be correlated :	
					db_cc[-1]['cpl'].append([ista1,ista2])
					ncpl+=1
		return db_cc, kwrong_dist 	


	def disp_path_repartition_per_file(self,dist,lvl=1) :#(ex,dist,npath_wrong_dist,lvl=1) :
		''' display information on the path that will be correlated'''
		# check that we have some paths. Otherwise exit
		#if lvl == 0 :
		#	return
		size_single_path = self.ex['size_single_path']
		if len(self.cc) == 0 : 
			dd.dispc('--------------','y','n')
			dd.dispc('  sorry no path were selected :/ Change your dist setting or set inner/outer to True','y','n')
			return 

		# print detailed number of path per h5 file 
		if lvl >=2 :
			dd.dispc('','c','b')
			dd.dispc('------------- paths repartition per h5 files details -------','c','b')
			for icc in range(0, len(self.cc)) :
				str_ = '  '+self.set.loc[icc,'name'].ljust(3)+' : '
				str_+= str(len(self.cc[icc])).rjust(5)+' couples of stations'
				#str_+= '  '+str(self.npath_wrong_dist)+' paths rej. bc of wrong dist '+str(dist)
				dd.dispc(str_,'c','n')


		# print the summary of the number of path per h5 file :
		if lvl >= 1 :
			dd.dispc('','c','b')
			dd.dispc('------------ paths repartition per h5 file summary ------------','c','b')
			dd.dispc('  we do not show pair of network with no common year  ','c','n')
			#k=-1
			for iset in range(0, len(self.cc)) :
				#k = k+ 1 
				croot = self.set.loc[iset,'name'][1:-5]
				if iset ==0 :
					npath = 0 #len(iset['cpl'])
					nfile = 0
					root_previous = croot 
					k=1
				if croot  == root_previous :
					npath = npath + len(self.cc[iset])
					nfile = nfile + 1 
				else : # on demarre un nouveau set ici :
					str_= '  '+root_previous.split('/')[-1].ljust(20)+' : '+str(npath).rjust(5)
					str_+='  paths splitted into '+str(nfile).rjust(5)+' h5 files'
					dd.dispc(str_,'c','n')
					root_previous = croot #self.set.loc[iset,'name'] 
					npath = len(self.cc[iset])
					nfile = 1 

			str_  = '  '+croot.split('/')[-1].ljust(20)
			str_ += ' : '+str(npath).rjust(5)+'  paths splitted into '
			str_ += str(nfile).rjust(5)+' h5 files'
			dd.dispc(str_,'c','n')

		# print the total	 number of path and the number of paths rejected because of wrong dist : 
		dd.dispc(' ','c','b')
		dd.dispc('------------ in summary : ------------','c','b')
		npath = 0 
		nfile = 0
		for iset in self.cc :
			nfile +=1 
			npath = npath + len(iset) 
		dd.dispc('  '+str(npath)+'  paths will be correlated into '+str(nfile)+' files','c','n')
		dd.dispc('  disk space : about ' +str(round(npath*size_single_path,2))+' Go (do not include metadata)','c','b')

#######################################################################################################################
#
#
#  DB STA
#   
########################################################################################################################

class DB_sta : 
	def __init__(self,dir_name,path_,date_) : 
		self.name={}
		self.name['h5']     = dir_name+'/db_sta.h5'
		self.name['lock']   = dir_name+'/db_sta.'+str(os.getpid())+'.lock'
		self.name['search'] = dir_name+'/db_sta.*'
		self.init = False                
		# db_sta already exist return without reading it. 
		if os.path.isfile(self.name['h5']) : 
			dd.dispc('  db_sta.h5 : do nothing','g','n')

		# construct it :
		if len(glob.glob(self.name['search'])) ==0 :  # is someone else working on db_sta 
			create_lock_file(self.name['lock'])       # no => create the lock file 
			self.construct_from_db_files(path_,date_)
			self.init = True                          # indicate that the object is loaded
			self.save()                               # save and remove lock file 


	def save(self):
		''' output the station info into h5 if it does not exist'''
		if os.path.isfile(self.name['h5']) :
			dd.dispc('  db_sta.h5 already exist => we do not need to save it','g','n')
			return
		# did this process we create a lock file with the right pid :
		if os.path.isfile(self.name['lock']) :
			dd.dispc('  we created the lock file => we create db_sta.h5','g','n')
			h5.save_class_data_to_h5(self,self.name['h5'])
			for iset in range(0,len(self.sta)) : 
				set_name = self.set.at[iset,'name']
				self.sta[iset].to_hdf(self.name['h5'],'/sta/'+set_name)
			# save the station file in a flat way to speed up the way we read CC : 
			sta_flat = pd.concat(self.sta,ignore_index=True)
			pyos.remove(self.name['lock'])


	def read(self) :
		''' import directly fe, cut_len, set, year, and add manually sta which is a list of dataframe'''
		h5.import_h5_data_to_class(self,self.name['h5']) 
		self.sta=[]
		for iset in self.set['name'] :
			self.sta.append(pd.read_hdf(self.name['h5'],'sta/'+iset))
		return 


	def construct_from_db_files(self,path_,date_) : 
		'''Create a dict of station sets and stations metadata + their availability year per year'''
		# create a list of years : 
		dd.dispc('  create '+self.name['lock']+' file and construct db_sta','g','n')
		year_list = []
		for idate_group in date_  :
			for idate in idate_group :
				cyear = UTCDateTime(idate).year 
				if not cyear in year_list :
					year_list.append(cyear)


		# construct the yearly availability frame
		set_name = []
		for iset in path_ : 
			set_name.append(iset.split('/')[-1]) 
		self.year = pd.DataFrame(columns=set_name,index=year_list)

		# initialize the data structure : 
		self.set = pd.DataFrame(columns=['name','path','year']) # one row per set 
		self.sta = []                                           # list of df ['name','id', 'lon', 'lat','elev','depth']
		set_= []                                                # temporary list of dict that will be converted to df
		# Now create a dict of set and stations with their per year availability : 
		nyear = len(self.year)
		# loop on set : 
		kset = -1
		for iset in path_  :
			kset = kset + 1 
			set_name = iset.split('/')[-1]
			set_.append({'path' : iset, 'name' : set_name}) # 'year_list' : np.zeros(nyear)})
			sta_list = []  # temporary list of dict that will be converted to df 
			sta_name = []  # temporary list of station name used to know if a station is already here or not
			kyear = -1 
			# loop on each year 
			for iyear in year_list :
				kyear = kyear + 1 
				db_name = iset+'/'+str(iyear)+'/db.pkl'
				# if the db file exist for this set/year open it and complete station list and set availability
				if os.path.isfile(db_name) :
					self.year.loc[iyear,set_name] = True
					#set_[-1]['year_list'][kyear]=True
					ff = open(db_name,'rb')
					db_iset = pickle.load(ff)
					ff.close()
					for ista in db_iset['sta'] :
						if not ista in sta_name :
							sta_name.append(ista)
							csta={}
							csta['lon']  = np.float64(db_iset['sta'][ista]['lon'])
							csta['lat']  = np.float64(db_iset['sta'][ista]['lat'])
							csta['elev'] = np.float32(db_iset['sta'][ista]['elev'])
							csta['depth']= np.float32(db_iset['sta'][ista]['depth'])
							csta['name'] = db_iset['sta'][ista]['name']
							csta['id']   = db_iset['sta'][ista]['kname'].replace('_','.')
							sta_list.append(csta)
			self.sta.append(pd.DataFrame(sta_list,columns=['name','id','lon','lat','elev','depth']))
		self.set = pd.DataFrame(set_,columns=['name','path']) #,'year_list'])
		self.fe      = db_iset['in_']['pp']['freq']
		self.cut_len = db_iset['in_']['pp']['cut_len'] # t1 en [s] 


	def disp_station_info(self) : 
		''' display station info on screen'''
		nset = self.set.shape[0]
		nsta =0
		for iset in range(0,nset) :
			nsta +=  self.sta[iset].shape[0] 
		ncpl = int((nsta*(nsta-1))/2 + nsta)
		dd.dispc('','c','n')
		dd.dispc('---------------stations info summary -----------------------','y','b')
		dd.dispc('  we got '+str(nset)+' sets of stations and a total of '+str(nsta)+' stations','y','n'); 	
		dd.dispc('  = '+str(ncpl)+' couples of stations if we correlate everything','y','n')


	def disp_full_station_list(self) :
		print('')
		dd.dispc('------------ full list of stations per -----------------','g','b')
		nset = self.set.shape[0]
		for iset in range(0,nset) :
			nsta = self.sta[iset].shape[0]
			dd.dispc('  '+self.set.at[iset,'name']+' :'+str(nsta).rjust(4)+' stations','g','n')
			pd.set_option('display.max_rows', self.sta[iset].shape[0]+1)
			print(self.sta[iset])

	def disp_station_list(self,lvl) :
		''' print for each set of stations the total number of stations found in db files '''
		print('')
		dd.dispc('------------ number of stations per set -----------------','g','b')
		nset = self.set.shape[0]
		for iset in range(0,nset) :
			nsta = self.sta[iset].shape[0]
			dd.dispc('  '+self.set.at[iset,'name']+' :'+str(nsta).rjust(4)+' stations','g','n')

	def is_there_a_common_year(self,kset1,kset2) :
		''' return True or False depending if there is a common year between the set #kset1 and #kset2'''
		nyear = len(self.year)
		for iyear in range (0,nyear) :
			if self.year.iloc[iyear,kset1] == True and self.year.iloc[iyear,kset2] == True :
				return True 
		return False


	def get_dist_and_azimuth(self,kset1,ksta1,kset2,ksta2)	:
		''' get distance and azimuth btw two stations from their set and station indice '''
		lon1 = self.sta[kset1].at[ksta1,'lon']
		lon2 = self.sta[kset2].at[ksta2,'lon']
		lat1 = self.sta[kset1].at[ksta1,'lat']
		lat2 = self.sta[kset2].at[ksta2,'lat']
		dist_,az,baz = gps2dist_azimuth(lat1,lon1,lat2,lon2)
		dist_ = dist_/1000 
		return dist_, az, baz,lon1,lat1,lon2,lat2 


	def _get_indice_first_station_in_each_set(self) :
		''' return a list containing the absolute indice of the 1st station of each set. Used by generate_db_file to convert relative indice to absolute one'''
		I_sta = [] 
		nsta = 0 
		for iset in range(0,len(self.sta)) :
			I_sta.append(nsta)
			nsta = nsta + len(self.sta[0])
		return I_sta 




##########################################################################################
##########################################################################################
#######    ###    ####         ####   ####   #######   ####         #####    #####    ####
#######    ###    ####         ####   ####   #######   ####         ######    ####    ####
#######    ###    #######   #######   ####   #######   #######   ##########   ###    #####
#######    ###    #######   #######   ####   #######   #######   ###########        ######
#######    ###    #######   #######   ####   #######   #######   ############      #######
#######    ###    #######   #######   ####   #######   #######   ############    #########
#######           #######   #######   ####        ##   #######   ###########    ##########
#######           #######   #######   ####        ##   #######   ##########    ###########
##########################################################################################
##########################################################################################

def create_lock_file(filename) : 
    ff=open(filename,'wb')
    pickle.dump(filename,ff)
    ff.close()


def	print_all_path (db_cc,db_sta) :
	for irow, iset in db_cc.set.iterrows() :
		iset1 = iset['I_set1']
		iset2 = iset['I_set2']
		filename = iset['name']
		for icpl, vcpl in db_cc.cc[irow].iterrows() : 
			sta_name1 = db_sta.sta[iset1]['id'][vcpl['I_sta1']]
			sta_name2 = db_sta.sta[iset2]['id'][vcpl['I_sta2']]

			print(filename+':   '+sta_name1+' - '+sta_name2)

