#general
import os , glob , h5py, pickle
#obspy
import obspy.core as obs 
from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth as gps2dist 

#scientific 
import numpy as np
import scipy.signal as signal 
from scipy.io import savemat as savemat 
import m_pycorr.mods.s2d as s2d 
import matplotlib.pyplot as plt
import matplotlib.patches as pt 
from matplotlib.collections import PatchCollection
#pycorr
import m_pycorr.mods.h5 as h5 
import m_pycorr.mods.dd as dd 
import m_pycorr.mods.lang as lang 
from m_pycorr_c3.live.c1_trace import c1_trace as c1_trace  

try : 
	import ipdb
except :
	pass 

class c3_live :
	def _get_filename_from_indice(self,kpath):
	 	# return the filename of the correlation #I 
		I_file = np.where(self.file_I[:,1] >= kpath)[0][0]
		return self.file[I_file]
	
	def c2_to_c1(self,kpath,kcmp,varargin_select_vs={},varargin={}) : 
		# selection par defaut des virtual source :
		in_vs= {} 
		in_vs['dist'] = [10, 1000] 
		in_vs['daz']  = 45 
		in_vs['rm_btw_sta'] = True 
		in_vs['single_source'] = False
		in_vs=lang.parse_options(in_vs,varargin_select_vs)
		#
		in_ = {}
		in_['savemat'] = False # only for testing purpose 
		in_['norm'] = True     # normalize the contribution of each vs b4 summing them ?
		in_['max_pos']=True    # always put the max in the pos side ?
		in_['tag'] = '' 
		in_ = lang.parse_options(in_,varargin)

		#determine directory name where the new C1 will be saved :
		out_dir = 'C3_to_C1__dist_'+str(in_vs['dist'][0])+'_'+str(in_vs['dist'][1])+'km'
		out_dir = out_dir +'__daz_'+str(in_vs['daz'])
		if in_['norm'] : 
			out_dir = out_dir +'__norm'
		if in_vs['rm_btw_sta'] :
			out_dir = out_dir +'__rm_btw_sta'
		out_dir = out_dir + in_['tag']

		# select virtual sources based on dist/az/...
		I1,I2, mat = self.select_virtual_source(kpath,kcmp,in_vs) 

		if in_vs['single_source'] : 
			if len(I1) >0 :
				I_min= np.argmin(mat['vs']['I1_daz'])
				I1 = np.array((I1[I_min],I1[I_min]))
			if len(I2) > 0 :
				I_min= np.argmin(mat['vs']['I2_daz'])					
				I2 = np.array((I2[I_min],I2[I_min]))

		# read the C3 :
		c3_vsmat_pp = self.read_c3_vsmat(kpath,kcmp,0) 
		c3_vsmat_nn = self.read_c3_vsmat(kpath,kcmp,1)

		# put the max always on the pos side :
		#if len(I1)>0 :
		if in_['max_pos'] :
			cc_pp_I1 = np.fliplr(c3_vsmat_pp.tr[I1,:])
			cc_nn_I1 = np.fliplr(c3_vsmat_nn.tr[I1,:])

		#if len(I2)>0 :
		cc_pp_I2 = c3_vsmat_pp.tr[I2,:]
		cc_nn_I2 = c3_vsmat_nn.tr[I2,:]

		# normalize each row before summing : 
		for itr in range(len(I1))  :
			max_ = max(abs(cc_pp_I1[itr,:])) 
			if max_ > 0 :
				cc_pp_I1[itr,:] = cc_pp_I1[itr,:] / max_ 
			max_ = max(abs(cc_nn_I1[itr,:]))
			if max_ > 0 :
				cc_nn_I1[itr,:] = cc_nn_I1[itr,:] / max_ 
		for itr in range(len(I2))  :
			max_ = max(abs(cc_pp_I2[itr,:]))
			if max_ > 0 :
				cc_pp_I2[itr,:] = cc_pp_I2[itr,:] / max_ 
			max_ = max(abs(cc_nn_I2[itr,:]))
			if max_ > 0 :
				cc_nn_I2[itr,:] = cc_nn_I2[itr,:] / max_ 

		# sum each matrix 
		cc_pp_I1 = cc_pp_I1.sum(axis=0)
		cc_pp_I2 = cc_pp_I2.sum(axis=0)
		cc_nn_I1 = cc_nn_I1.sum(axis=0)
		cc_nn_I2 = cc_nn_I2.sum(axis=0)

		cc_sum = cc_pp_I1 + cc_pp_I2 + cc_nn_I1 + cc_nn_I2

		# complete the .mat ict containing all info 
		mat['cc_pp_I1'] = cc_pp_I1
		mat['cc_pp_I2'] = cc_pp_I2
		mat['cc_nn_I1'] = cc_nn_I1
		mat['cc_nn_I2'] = cc_nn_I2
		mat['cc_sum']   = cc_sum
		mat['time'] = self.t
		mat['out_dir'] = out_dir           # to know where to save plots :-)

		# save is requested : 
		if in_['savemat'] :
			try :
				c1 = self.read_c1_single_pair(kpath,kcmp)
				mat['c1'] = c1.tr
				mat['c1_time'] = c1.time
			except :
				pass 
			out_dir_mat = 'test/'+out_dir 
			try :
				os.makedirs(out_dir_mat) 
			except :
				pass 
			savemat(out_dir_mat+'/'+mat['name'][0:-3].replace('.','_')+'.mat',mat)
		return mat 


	def select_virtual_source(self,kpath,kcmp,varargin={}) : 
		# renvoit les indices des virtual source du cote sta1, sta2, et l'ensemble 

		#config type 1:
		#in_= {} 
		#in_['dist'] = [10, 1000] 
		#in_['daz']  = 45 
		#in_['norm'] = False
		#in_['rm_btw_sta'] = True  
		#in_=lang.parse_options(in_,varargin)

		#config type 2:
		in_= {} 
		in_['dist'] = [10, 1000] 
		in_['daz']  = 360
		in_['rm_btw_sta'] = False
		in_=lang.parse_options(in_,varargin)

		# prepare the vs dict containing informations about virtual sources : 
		vs   = {}
		nvs  = len(self.vs['lat'])
		vs['dist1'] = np.zeros((nvs)) # distance to sta1 
		vs['dist2'] = np.zeros((nvs)) # distance to sta2 
		vs['distc'] = np.zeros((nvs)) # distance to center 
		vs['az1']   = np.zeros((nvs)) # az vs to sta1 
		vs['az2']   = np.zeros((nvs)) # az to sta2 
		vs['azc']  = np.zeros((nvs))  # az from vs to middle of sta1,sta2
		vs['bazc']  = np.zeros((nvs)) # az  from vs to middle of sta1,sta2.

		# for convenience : extract the lon, lat of the two stations 
		lon = self.lon[kpath]
		lat = self.lat[kpath]
		lon_mean = np.mean(lon)
		lat_mean = np.mean(lat) 

		# compute distanc/az/baz of each Vs with respect to the sta1, sta2 and the center of sta1,sta2 :
		for ivs in range(nvs) : 
			[vs['dist1'][ivs], vs['az1'][ivs], baz] = gps2dist(self.vs['lat'][ivs],self.vs['lon'][ivs],lat[0],lon[0])
			[vs['dist2'][ivs], vs['az2'][ivs], baz] = gps2dist(self.vs['lat'][ivs],self.vs['lon'][ivs],lat[1],lon[1])
			[vs['distc'][ivs], vs['azc'][ivs], vs['bazc'][ivs]] = gps2dist(self.vs['lat'][ivs],self.vs['lon'][ivs],lat_mean,lon_mean)

		vs['dist1'] = vs['dist1']/1000 
		vs['dist2'] = vs['dist2']/1000 
		vs['distc'] = vs['distc']/1000 

		# indice of virtual sources closer to sta1 or sta2 :
		I_sta1 = np.where(vs['dist1'] < vs['dist2'])[0]   # indice of virtual sources closer to sta1 
		I_sta2 = np.where(vs['dist2'] <= vs['dist1'])[0]  # indice of virtual sources closer to sta2


		# list virtual source at more than 90 degree from nearest station i.e not between the stations : 
		az  = self.az[kpath] +90 
		baz = self.baz[kpath] +90 
		if baz > 360 :
			baz = baz-360 
		if az > 360 : 
			az = az - 360  
		I1=np.where((vs['az1'] < az-10) & (vs['az2'] < az-10))[0]
		I2=np.where((vs['az1'] > az+10) & (vs['az2'] > az+10))[0]
		I3=np.where((vs['az1'] < baz-10) & (vs['az2'] < baz-10))[0]
		I4=np.where((vs['az1'] > baz+10) & (vs['az2'] > baz+10))[0]
		I_not_btw = np.intersect1d(np.concatenate((I3,I4)),np.concatenate((I1,I2)))

		# list virtual source at the right distance from sta1, sta2 and the center of sta1,sta2 
		I1 = np.where((vs['dist1'] >= in_['dist'][0]) & (vs['dist1'] <= in_['dist'][1]))[0]
		I2 = np.where((vs['dist2'] >= in_['dist'][0]) & (vs['dist2'] <= in_['dist'][1]))[0]
		I3 = np.where((vs['distc'] >= self.dist[kpath]/2))
		I_dist = np.intersect1d(I1,I2);
		I_dist = np.intersect1d(I_dist,I3);

		# list virtual sources so vs->sta1 or vs-> sta2 has an acceptable azimuth : 
		I0 = np.where(abs(vs['az1'][I_sta1] - self.az[kpath]+360) <= in_['daz'])[0]
		I1 = np.where(abs(vs['az1'][I_sta1] - self.az[kpath]) <= in_['daz'])[0]
		I2 = np.where(abs(vs['az2'][I_sta2] - self.baz[kpath]) <= in_['daz'])[0]
		I3 = np.where(abs(vs['az2'][I_sta2] - self.baz[kpath]+360) <= in_['daz'])[0]
		I_az = np.concatenate((I_sta1[I0],I_sta1[I1],I_sta2[I2],I_sta2[I3]))
		#I_close_sta1=I1 
		#I_close_sta2=I2
		# apply all criterias => I_sel = indice of Vs selected  
		I_sel = np.intersect1d(I_az,I_dist)
		if in_['rm_btw_sta'] == True :
			I_sel = np.intersect1d(I_sel,I_not_btw)

		I1 = np.intersect1d(I_sel,I_sta1) # indice of selected Vs closer to sta1 
		I2 = np.intersect1d(I_sel,I_sta2) # --------------------------------sta2 

		

		# create a mat dict that can be saved into a .mat file for plotting
		mat = {}
		mat={}
		mat['vs'] = self.vs 
		mat['vs']['I_close_sta1']=I_sta1
		mat['vs']['I_close_sta2']=I_sta2
		mat['vs_info']=vs 
		mat['vs']['I_az'] = I_az           # indice of vs having the right az 
		mat['vs']['Idist'] = I_dist        # indice of vs at the right distance
		mat['vs']['I_not_btw'] = I_not_btw # indice of vs not btw the 2 stations
		mat['vs']['I_sel']     = I_sel     # indice of vs selected
		mat['name'] = 'dist_'+str(np.round(self.dist[kpath]))+'km_'+str(self.id[kpath][0],'utf8')+'_'+str(self.id[kpath][1],'utf8')+'.mat'
		mat['lon'] = self.lon[kpath]       # lon, lat, dist of the station couple 
		mat['lat'] = self.lat[kpath]       #
		mat['dist'] = self.dist[kpath]     #
		mat['az'] = self.az[kpath]
		mat['baz'] = self.baz[kpath]
		mat['vs']['I1'] = I1               # indice of Vs selected closer to sta 1 
		mat['vs']['I2'] = I2               # indice of Vs selected closer to sta 2
		mat['vs']['I1_daz'] = abs(vs['az1'][I1] - self.az[kpath])
		mat['vs']['I2_daz'] = abs(vs['az2'][I2] - self.baz[kpath])
		return I1, I2, mat
		

	def read_c3_vsmat(self,kpath,kcmp,kc3_type) :
		# return a c1_trace object + plain matrix  + vs position ?
		dset_name ='/c3_mat/'+str(self.id[kpath,0],'utf8')+'/'+str(self.id[kpath,1],'utf8')+'/'+str(self.in_['cmp'][kcmp],'utf8')
		dset_name ='/c3_mat/'+str(self.id[kpath,0],'utf8')+'/'+str(self.id[kpath,1],'utf8')+'/'+str(self.in_['cmp'][kcmp],'utf8')
		dset_name = dset_name+'/'+str(self.in_['c3_type'][kc3_type],'utf8')
		filename = self._get_filename_from_indice(kpath) 
		tr = self.read_single_pair(filename,dset_name,kpath,kcmp)
		return tr

	def read_c1_single_pair(self,kpath,kcmp) :
		# always zz 
		dset_name ='/c1/'+self.id[kpath,0]+'/'+self.id[kpath,1]+'/ZZ'
		filename = self._get_filename_from_indice(kpath) 
		tr = self.read_single_pair(filename,dset_name,kpath,kcmp)
		if len(tr.tr > 0) : # if we found a c1 :
			tr.time = self.t_c1a
		return tr

	def read_c3_single_pair(self,kpath,kcmp,kc3_type) : #always work on /c3 
		dset_name ='/c3/'+self.id[kpath,0]+'/'+self.id[kpath,1]+'/'+self.in_['cmp'][kcmp]
		dset_name = dset_name+'/'+self.in_['c3_type'][kc3_type]
		filename = self._get_filename_from_indice(kpath) 
		return self.read_single_pair(filename,dset_name,kpath,kcmp)

	def read_single_pair(self,filename,dset_name,kpath,kcmp) :
		# read the trace dset_name in filename and return c1_trace object
		# if trace does not exist tr['tr'] is empty
		# return c3 time vector even if we attempt to read a c1 (corrected in read_c1_single_pair)
		ff = h5py.File(filename,'r') 
		tr = {}
		if dset_name in ff :
			tr['tr']=ff[dset_name][:]
		else 	:
			tr['tr'] = [] 
		if len(self.md_c['date1']) > 0 :
			tr['date1'] = self.md_c['date1'][0]
			tr['date2'] = self.md_c['date2'][-1]
		tr['title'] = dset_name
		tr['az']    = self.az[kpath] 
		tr['baz']   = self.baz[kpath]
		tr['cmp']   = self.in_['cmp'][kcmp]
		tr['lon']   = self.lon[kpath]
		tr['lat']   = self.lat[kpath]
		tr['id']    = self.id[kpath]
		tr['dist']  = self.dist[kpath] 
		tr['elev']  = self.elev[kpath]
		tr['file']  = filename
		tr['time']  = self.md_c['t']
		tr['tau']   = self.md_c['tau']
		return c1_trace(tr)


	def __init__(self,dir_name) :
		db_name = dir_name+'/'+'db.pkl' 
		if not os.path.isfile(db_name) :
			dd.dispc('  '+db_name+' doest not exist, lets create it','c','b')
			make_db(dir_name,db_name)
		
		dd.dispc('  loading '+db_name,'c','b')
		db=load_pkl(db_name)
		for ikey in db.keys() : 
			setattr(self,ikey,db[ikey])



#------------------------------------------------------
#   function that are outside the class 
#------------------------------------------------------
def make_db(dir_name,db_name) :
	# get xcorr_*h5 file list :
	file_list = glob.glob(dir_name+'/xcorr*h5');
	file_list.sort() 
	nfile=len(file_list) 
	# init structure : 
	db={}
	db['lat'] = [] 
	db['lon'] = [] 
	db['elev']= [] 
	db['id']  = []
	db['file_I']=np.zeros((nfile,2))
	
	# 1. loop on all xcorr*h5 file and read their metadata
	npath_prev = 0 
	k=0
	for ifile in file_list : 
		ff = h5py.File(ifile,'r')
		db['lat'].extend(ff['md']['lat'][:].tolist())
		db['lon'].extend(ff['md']['lon'][:].tolist())
		db['elev'].extend(ff['md']['elev'][:].tolist())
		db['id'].extend(ff['md']['id'][:].tolist())
		db['file_I'][k,:] = [npath_prev, len(db['lat'])-1]
		npath_prev = len(db['lat'])
		ff.close()
		k=k+1

	# converting list into array : 
	db['id']  = np.array(db['id'])
	db['lon'] = np.array(db['lon']) 
	db['lat'] = np.array(db['lat'])
	db['elev']= np.array(db['elev'])
		
	#computing distance and baz btw each station pair : 
	npath = len(db['lat']) 
	db['dist']= np.zeros((npath,1))
	db['az']  = np.zeros((npath,1))
	db['baz'] = np.zeros((npath,1))
	for ipath in range(0, npath) :
		[dist,az,baz] = gps2dist(db['lat'][ipath][0],db['lon'][ipath][0],db['lat'][ipath][1],db['lon'][ipath][1])
		db['dist'][ipath] = dist/1000. 
		db['az'][ipath]   = az 
		db['baz'][ipath]  = baz 
	# adding md_c, vs , in :
	ff = h5py.File(file_list[0])
	db['md_c']= h5.read_group_as_dict(ff['md_c'])
	db['in_'] = h5.read_group_as_dict(ff['in_']) 
	vs = h5.read_group_as_dict(ff['vs']) 
	nvs = len(vs['lat'])
	#for ivs in range (nvs) :
	#	vs['lat'][ivs] = vs['lat'][ivs][0]
	#ipdb.set_trace()
	# we need to make it flat :/ 
	db['vs'] = {}
	for ikey in ['lon','lat','elev'] :#vs.keys():
		db['vs'][ikey] = np.zeros((nvs))#
		print(ikey)
		for ivs in range (nvs)  : 
			db['vs'][ikey][ivs] = vs[ikey][ivs][0]
	db['vs']['id'] = [] 
	for ivs in range(nvs) : 
		db['vs']['id'].append(vs['id'][ivs][0])

	# create some shortcut to mimic the matlab structure :
	db['t'] = db['md_c']['t']
	db['t_c1a'] = db['md_c']['t_c1a']
	db['t_c1b'] = db['md_c']['t_c1b']
	db['tau']  = db['md_c']['tau']
	db['date1'] = db['md_c']['date1']
	db['date2'] = db['md_c']['date2']
	db['file']  = file_list 
	save_as_pkl(db_name,db)


def load_pkl(filename) : 
	#dd.dispc('  loading '+filename,'y','b')
	print(filename)
	ff=open(filename,'rb') 
	db=pickle.load(ff)
	ff.close() 
	return db 

def save_as_pkl(filename,db) :
	ff=open(filename,'wb')
	pickle.dump(db,ff)
	ff.close()

