
import m_pycorr.analyse.mods.cc_h5 as cc 
import pycorr_mods.dd as dd 
import pycorr_mods.lang as lang 
import h5py 
import numpy as np 
import os,glob, pickle, random 
from obspy.geodetics.base import gps2dist_azimuth as get_distance

try : 
	import ipdb 
except	:
	pass 

# TODO : 
# -ajouter dist,az,baz, cmp , ref_nstack, cc_nstack
# -ajouter file file_I 
# - temporairement : try/except lorsqu'on concatene ref_nstack. Cela provient 
#   du faite que le code de rotation de Pierre ne remplit pas ce champ 
# 	
# HYPOTHESES :
# -- toutes les metadonnes des correlations [/md_c/lon ...] sont de taille [npath,2]
# -- pas plus de un million de group et un milliard de trajets/groupe : 
# -- un ls pendant que le code tourne doit renvoyer les groupes/paths dans le meme ordre que dans les fichiers 
#    de correlations 
# -- dans les fichiers de correlations : tous les trajets doivent avoir des donnees de la meme taille i.e 
#    des correlations horaires/journalieres/reference de taille identique  
#
# FORMALISME RETOUR DES INFOS : les fonction processant les donnees renvoient : 
# -dvv  = dictionnaire du meme nom que la fonction e.g dvv, dt, dispersin,... 
# -in_  = dictionnaire contenant parametre input de la fonction appelle (dvv)
# -md_c = metadonnees communes a tous les trajets testes (pour dvv les data, pour dispersion les periods,...) 
#
# DECOMPOSITION DE LA BOUCLE ET DE LA DISTRIBUTION DES PROCESSUS :
# C1/pydb_dvv_XX.h5 [.pid.lock .pid.h5]                                       => pydb_name
# C1/pydb_dvv_XX/xcorr_set1_set2  [.pid.lock .pid.h5 /]                       => pydb_xcorr 
# C1/pydb_dvv_XX/xcorr_set1_set2/group_00001  [.pid.lock .pid.h5 /]           => pydb_group /
# C1/pydb_dvv_XX/xcorr_set1_set2/group_00001/path_001_005 [.pid.lock .pid.h5] => pydb_path
#
# structure fichier hdf5 :
# /md       : metadonnees identique aux correlation 
# /md_c     : metadonnes communes aux correlations 
# /dvv/in_  : parametres input des mesures de dv/v 
# /dvv/md_c : metadonnees communes aux dv/v 
# /dvv/corrcoeff/sta1/sta2/ZZ 
# /dvv/dvv/sta1/sta2/ZZ
#
# c1 => everything that relates to correlations stored in h5 format 

def main_loop(in_dir,fh,in_fh,in_ml) :
	''' in_dir : directory where correlations are stored i.e 'C1_01_xcorr___coher'
		fh     : pointer to a function processing correlations computed between a single pair of station
		in_fh  : input parameters for fh [dict]
		in_ml  : input parameters for main_loop_h5
	''' 
	in_={}
	in_['process_1st_group_only']= False  #to be implemented cut the group number in define_group ?
	in_['load_ref']       = True  
	in_['load_cc']        = True  
	in_['npath_per_core'] = 10 
	in_['npath_per_group']= 5000
	in_['cmp']             =['ZZ']
	in_=lang.parse_options(in_,in_ml)
	dd.dd(in_)

	#define the final filename : C1/pydb_dvv_XX. [h5 .pid.lock .pid.h5]
	pydb_names = ex_define_name(in_dir+'/'+in_fh['tag']) 
	dd.dispc(pydb_names['root'],'y','b')
	if ex_has_everything_been_done(pydb_names,2) : return 

	#list C1/xcorr_set1_set2.h5 files and randomize their order :
	c1_list = c1_list_correlation_files(in_dir)
	c1_list_random = c1_list_correlation_files(in_dir)
	if not in_['process_1st_group_only'] :
		random.shuffle(c1_list_random)

	# define C1/pydb_dvv/xcorr_set1_set2 root name 
	xcorr_root= ex_define_xcorr_root(pydb_names['root'],c1_list_random)
	# loop on C1/xcorr_set1_set2.h5 files :
	for icc in c1_list_random :
		pydb_xcorr = ex_define_name(xcorr_root[icc])
		if ex_has_everything_been_done(pydb_xcorr) : continue 
		process_this_xcorr_file(icc,pydb_xcorr,fh,in_fh,in_)
	if ex_should_we_concatenate(pydb_names,len(c1_list_random)) :
		concatenate(pydb_names,0,prefix='/'+fh.__name__+'/')
		add_metadata_from_xcorr_file(pydb_names['h5'],c1_list)


#---------------------------------------------
def process_this_xcorr_file(c1_name,pydb_xcorr,fh,in_fh,in_):
	mkdir(pydb_xcorr['root'])
	#get a dictionnary groups name and indice of path they will contain : 
	groups =define_group(c1_name,pydb_xcorr['root'],in_['npath_per_group'],in_['process_1st_group_only']) 
	ngroup=len(groups)

	if in_['process_1st_group_only'] == True : 
		#group_key = sorted(groups.iteritems())
		group_key = sorted(iter(groups.items()))
	else :	
		group_key = iter(groups.items())
		#group_key = groups.iteritems()

	for igroup_name, igroup in group_key : 
		#get the output file name of this group :/c1/pydb_dvv/xcorr_set1_set2/group_000 [.h5...]
		pydb_group = ex_define_name(igroup['root'])
		if ex_has_everything_been_done(pydb_group,8) : continue #pas de .h5, .pkl, .lock
		mkdir(pydb_group['root'])
		path_list=define_path_list(igroup,in_['npath_per_core'])
		process_this_group(igroup,path_list,fh,in_fh,in_)
		if  ex_should_we_concatenate(pydb_group,len(path_list)) : 
			concatenate(pydb_group,8)
	if  ex_should_we_concatenate(pydb_xcorr,ngroup) : 
		concatenate(pydb_xcorr,4)


def process_this_group(igroup, path_list,fh,in_fh,in_):
	dd.dispc('    '+igroup['root'],'c','b')
	mkdir(igroup['root'])
	#loop on path in a random order : 	
	path_list_random=path_list.values() 
	if in_['process_1st_group_only'] :
		random.shuffle(path_list_random)
	for ipath in path_list_random  : 
		pydb_path = ex_define_name(ipath['root'])
		if ex_has_everything_been_done(pydb_path,12) : continue 
		create_lock_file(pydb_path['lock_pid'])
		dd.dispc('        '+pydb_path['root'],'c','n')
		result=process_these_paths(igroup['c1_name'],ipath,fh,in_fh,in_)
		save_as_hdf5(pydb_path,result,in_['cmp'])
		remove(pydb_path['lock_pid'])


def process_these_paths(c1_name,kpath,fh,in_fh,in_) : 
	ff=h5py.File(c1_name,'r')
	cmp_list=in_['cmp'] #ff['/md_c/cmp'][:]
	result={}
	result['h5_path']=[]
	for icmp in cmp_list : 
		result[icmp]=[]
		for ipath in range(kpath['I1'],kpath['I2']+1,1) :
			id_=ff['/md/id'][ipath]
			path_in_h5='/'+id_[0].decode()+'/'+id_[1].decode()+'/'+icmp
			dd.dispc(path_in_h5,'c','n')
			data={}
			if in_['load_ref'] : 
				try :
					data['ref']    =ff['/ref'+path_in_h5][:]
				except	:
					ipdb.set_trace()
			if in_['load_cc'] :
				data['cc_mat'] =ff['/cc'+path_in_h5]   #PROTECT
			data['t']    =ff['/md_c/t'][:]
			data['date1']=ff['/md_c/date1_ordinal'][:]
			data['date2']=ff['/md_c/date2_ordinal'][:]
			data['lon']=ff['/md/lon'][ipath][0:2]
			data['lat']=ff['/md/lat'][ipath][0:2]
			lon0 = data['lon'][0]
			lon1 = data['lon'][1]
			lat0 = data['lat'][0]
			lat1 = data['lat'][1]
			data['dist']=get_distance(lat0,lon0,lat1,lon1)[0]/1000
			data['id'] = ff['/md/id'][ipath]
			result_for_this_path, md_c, in_fh_used =fh(data,in_fh)
			result[icmp].append(result_for_this_path)
			if icmp == cmp_list[0]:
				result['h5_path'].append('/'+ff['/md/id'][ipath][0].decode()+'/'+ff['/md/id'][ipath][1].decode())
	ff.close()
	result['md_c']  = md_c
	result['in_']   = in_fh_used 
	result['in_ml'] = in_  
	return result



#---------------------------------------------------------------
#
#                 PROCESS THIS *  SUBFUNCTIONS 
#
#--------------------------------------------------------------

def define_group(c1_name,xcorr_root,npath_per_group,process_1st_group_only) : 
	''' called by process_this_xcorr_file. Define 
		-name of the C1/pydb/dvv/xcorr_set1_set2/group_* directories where
		-Indice of the 1st and last station pair belonging to each group
	'''
	dd.dispc(c1_name,'w','b')
	ff=h5py.File(c1_name,'r') 
	npath = len(ff['/md/lon'])
	gr={}
	k=-1
	for ipath in range(0,npath,npath_per_group) :
		k=k+1
		gr_name='group_'+str(k).zfill(6)
		gr_path=xcorr_root+'/group_'+str(k).zfill(6)
		gr[gr_name]={}
		gr[gr_name]['I1']=ipath 
		gr[gr_name]['I2']=min(ipath+npath_per_group-1,npath-1) 
		gr[gr_name]['c1_name']=c1_name
		gr[gr_name]['root']=gr_path
	ff.close()

	#randomize group keys : 
	if not process_1st_group_only :
		print('xx')
		gr_keys =list(gr.keys())
		random.shuffle(gr_keys)
		gr_rand={}
		for ikey in gr_keys : 
			gr_rand[ikey]=gr[ikey]
	else :
		gr_rand = gr 
	return gr_rand


def define_path_list(igroup,npath_per_core) :
	''' called by process this group : 
		-for a given group C1/pydb_dvv/xcorr_set1_set2/group_0000/ define 
		the name of each path file and the indice of 1st and last path they will 
		contain
	'''
	path={}
	for ipath in range(igroup['I1'],igroup['I2']+npath_per_core,npath_per_core) : 
		I1=ipath 
		I2=min(ipath+npath_per_core-1,igroup['I2'])
		if I2 < I1 : continue #je craque :
		path_name=str(I1).zfill(9)+'_'+str(I2).zfill(9)
		path_out=igroup['root']+'/path_'+str(I1).zfill(9)+'_'+str(I2).zfill(9)
		path[path_name]={}
		path[path_name]['root']=path_out
		path[path_name]['I1']=I1
		path[path_name]['I2']=I2
	return path 



#-------------------------------------------------------------------------
#
#                                HDF5 fUNCTIONS 
#
#------------------------------------------------------------------------------

def add_metadata_from_xcorr_file(pydb_file, c1_list) :
	''' pydb_file = C1_xcorr/pydb_dvv.h5 
		c1_list   = C1_xcorr/xcorr_set1_set2.h5 file list
		we copy and concatenate all c1_list files metadata ['/md/lon,lat,id,...']
		into the C1/xcorr/pydb_dvv.h5 file
	'''
	print('xx')
	ff=h5py.File(pydb_file,'a')
	#copy des metadata common to all correlations :  
	fx=h5py.File(c1_list[0],'r')
	ff.copy(fx['/md_c'],'/md_c') 
	ff.copy(fx['/in_'],'in_')

	#get information on each kind of metadata ([/md/lon, lat,...]) : name + dtype 
	md_list = [md_list for md_list in fx['/md'].keys()]
	#	md_list = fx['/md'].keys()
	md_dtype={}
	for imd in md_list :
		md_dtype[imd]=fx['/md/'+imd].dtype
	fx.close()
	#ici il y'a un bug possible : on se sert du premier fichier xcorr pour evaluer le dataype de md.id 
	#or il peut arriver dans le 1er fichier les id fasse par exemple 10 caractere, alors que dans les 
	#fichiers suivant on puisse en trouver de plus grand. 
	#=>on se met une marge de 2 caractere :
	ndtype=str(int(str(md_dtype['id'])[2:])+2)
	#print(ndtype)
	md_dtype['id']=np.dtype('S'+ndtype)

	#we need to know how many path we have to initialize the concatenated metadata :
	npath=0
	for ic1 in c1_list : 
		fx=h5py.File(ic1,'r')
		npath=npath + len(fx['/md/lon'])
		fx.close()

	#initialize metadata (assuming there are all 2D)

	try : # for C2 converted to C1 :/
		ncmp_cc = len(ff['/md_c']['cmp'])
	except	:
		ncmp_cc = 1 
	nfile   = len(c1_list)
	ff.create_dataset('/ref_nstack',(ncmp_cc,npath))
	ff.create_dataset('/file_I',(nfile,2))
	ff.create_dataset('/file',data=np.array(c1_list,dtype='S'))
	for imd in md_list : 
		ff.create_dataset('/md/'+imd,(npath,2),dtype=md_dtype[imd])
		
	#copy metadata from xcorr files to the pydb file : 
	I1=0
	kc1=0
	for ic1 in c1_list : 
		fx=h5py.File(ic1,'r')
		#print ic1 
		#print len(fx['/md/lon'])
		I2 = I1 + len(fx['/md/lon'])
		#print I1,I2
		ff['/file_I'][kc1,:]=[I1,I2-1]
		try :
			ff['/ref_nstack'][:,I1:I2]=fx['/ref_nstack']
		except :
			pass 
		for imd in md_list :
			ff['/md/'+imd][I1:I2,:]=fx['/md/'+imd]
		I1=I2
		kc1=kc1+1
		fx.close()
	# -ajouter dist,az,baz, cmp , ref_nstack
	# -ajouter file file_I 
	fmd=ff['/md']
	fdist= fmd.create_dataset('dist',(npath,1))
	faz  = fmd.create_dataset('az',(npath,1))
	fbaz = fmd.create_dataset('baz',(npath,1))
	for ipath in range(0,npath,1) :
		lat0 = fmd['lat'][ipath][0]
		lat1 = fmd['lat'][ipath][1]
		lon0 = fmd['lon'][ipath][0]
		lon1 = fmd['lon'][ipath][1]
		fdist[ipath],faz[ipath],fbaz[ipath]=get_distance(lat0,lon0,lat1,lon1) 
	fdist[:]=fdist[:]/1000.

	ff.close()

def concatenate(filenames,tab=4,prefix='') :
	''' concatenate a set of temporary h5 file located in the filenames['root'] dir into a new h5 file
	   (e.g C1/pydb_dvv/xcorr_set1_set2/*h5 => C1/pydb_dvv/xcorr_set1_set2.h5)
	    and remove all previous h5 files 
	'''
	create_lock_file(filenames['lock_pid'])
	dd.dispc(''.ljust(tab)+'concatenating '+filenames['root'],'y','n')
	ff=h5py.File(filenames['h5_pid'],'w')
	
	h5_list = glob.glob(filenames['root']+'/*.h5')
	h5_list.sort()
	for ih5 in h5_list :
		fh5 = h5py.File(ih5,'r')
		key_list = [key_list for key_list in fh5.keys()]
		#key_list=fh5.keys()
		key_list.remove('in_')
		key_list.remove('md_c')
		key_list.remove('in_ml')
		if ih5 == h5_list[0] :  # 1ere iteration :
			ff.copy(fh5['in_'],prefix+'in_')
			ff.copy(fh5['in_ml'],prefix+'in_ml')
			ff.copy(fh5['md_c'],prefix+'md_c')
		for ikey in key_list : 
			for ista1 in fh5[ikey].keys():
				for ista2 in fh5[ikey][ista1].keys():
					try : # can crash if the same path appears twice (when merging different set of CC)
						ff.copy(fh5[ikey][ista1][ista2],prefix+ikey+'/'+ista1+'/'+ista2+'/')  
					except	:
						pass
		fh5.close()
	ff.close()
	rename(filenames['h5_pid'],filenames['h5'])
	remove_dir(filenames['root'])
	remove(filenames['lock_pid'])


def save_as_hdf5(filename,data,cmp_list) : 
	''' save measurements corresponding to a set of path (C1/pydb_dvv/xcorr_set1_set2/group_xx/path_xx)
		into a hdf5 file having a structure :
		/md_c : metadata common to all measurements i.e independant of the station pair 
		/in_  : input parameters used by the function that processed the correlations 
		/param1/sta1/sta2/ZZ : result of the measurements 
		/param2/sta1/sts2/ZZ : -- 
	'''	
	ff=h5py.File(filename['h5_pid'],'w')
	h5_save_dict_keys_as_group(ff,{'md_c' :data['md_c'], 'in_' : data['in_'], 'in_ml' : data['in_ml']})
	key_list = data[cmp_list[0]][0].keys()
	npath = len(data[cmp_list[0]])
	for icmp in cmp_list : 
		for ipath in range(0,npath) :
			for ikey in key_list : 
				h5_path='/'+ikey+data['h5_path'][ipath]+'/'+icmp
				ff.create_dataset(h5_path,data=data[icmp][ipath][ikey])
	ff.close()
	rename(filename['h5_pid'],filename['h5'])


def h5_save_dict_keys_as_group(h5,dict_,prefix=''):
    for ikey in dict_ :
        h5.create_group(prefix+'/'+ikey)
        for dname, dset in iter(dict_[ikey].items()):
            if type(dset)==list and len(dset) > 0 and type(dset[0]) == dict :
                    h5_copy_dict_to_group({dname : dset[0]},h5,prefix='/'+ikey)
            else :
                try :
                    h5.create_dataset(prefix+'/'+ikey+'/'+dname,data=dset)
                except :
                	pass
    return h5 


#-------------------------------------------------------------------------
#
#           EX FUNCTIONS CONTROLLING THE FLOW OF THE CODE
#
#--------------------------------------------------------------------------

def ex_should_we_concatenate(filenames,nfile) :
	#c'est pas entrain d'etre concatene : =has_everything_been done 
	if ex_has_everything_been_done(filenames) : return False 
	#il ya le bon nombre de fichiers.h5  (sachant qu'il peut y avoir des fichiers temporaires ds la list .pid.h5)
	if len(glob.glob(filenames['root']+'/*.h5')) != nfile : return False
	#il ya pas de fichier /.lock
	if len(glob.glob(filenames['root']+'/*.lock')) > 0: return False
	return True 


def ex_define_xcorr_root(root,c1_list) : 
	pydb_xcorr={}
	for icc in c1_list :
		pydb_xcorr[icc]=root+'/'+icc[0:-3].split('/')[1]
	return pydb_xcorr


def ex_has_everything_been_done(filenames,tab=4) : 
	tab=''.ljust(tab)

	if os.path.isfile(filenames['h5']) : 
		dd.dispc(tab+filenames['h5']+' already exist','r','n')
		return True

	if len(glob.glob(filenames['search'])) >0 :
		dd.dispc(tab+filenames['search']+' already exist','r','n')
		return True
	return False 


def ex_define_name(root) : 
	pid=str(os.getpid())
	name={}
	name['root']    =root 
	name['h5']      =name['root']+'.h5' 
	name['h5_pid']  =name['root']+'.'+pid+'.h5'
	name['lock_pid']=name['root']+'.'+pid+'.lock'
	name['search']  =name['root']+'.*'
	return name 




#---------------------------------------------
#
#             MISC FUNCTIONS
#
#----------------------------------------------------

def c1_list_correlation_files(in_dir) : 
	c1_list=glob.glob(in_dir+'/[xv]*.h5')
	c1_list.sort()
	return c1_list 


def create_lock_file(filename) :
	ff=open(filename,'wb')
	pickle.dump(filename,ff)
	ff.close()



#----------------------------------------------------------------------------
#
#                           OS FUNCTION WITH PROTECTION 
#
#-----------------------------------------------------------------------------

def mkdir(dir_name) : 
	if not os.path.isdir(dir_name) : 
		try :
			os.makedirs(dir_name)
		except :
			pass

#---------------------
def remove(filename) : 
	try :
		os.remove(filename)
	except :
		pass 

#----------------------------
def remove_dir(dir_name) : 
	file_list=glob.glob(dir_name+'/*')
	for ifile in file_list  :
		try :
			os.remove(ifile)
		except :
			pass 
	try :
		os.removedirs(dir_name)
	except :
		pass 

def rename(src,dest) :
	try : 
		os.rename(src,dest)
	except 	:
		pass 










