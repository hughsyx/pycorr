################################################
# h5.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


import pandas as pd 
import h5py 
import ipdb


# WORK ON A CLOSE H5 FILE + A CLASS INSTANCE ON WHICH WE READ/WRITE :
def save_class_data_to_h5(self,filename) : 
	''' save dataframe, dict and scalar of the class instance "self" to a h5file. 
		We cannot each element of list of dict/dataframe to a different group since we have no way to name the group'''
	item = self.__dict__ .keys()
	ff = h5py.File(filename,'a') 
	for iname in item : 
		data = self.__dict__[iname]
		if type(data) == pd.DataFrame :
			data.to_hdf(filename,iname)
		elif type(data) == dict : 
			copy_dict_to_group(ff,{iname : data},'')			
		elif type(data) == list  :
			if len(data) == 0 :
				continue
		elif type(data) == float or type(data)==int : 
			ff.create_dataset(iname,data=data) 
	ff.close()


def import_h5_data_to_class(self,filename) : 
	ff =h5py.File(filename,'r') 
	for ikey in ff.keys() :
		# on a un group 
		if type(ff[ikey]) == h5py._hl.group.Group :
			# est-ce que c'est un dataframe ecrit par pandas : 
			if 'axis0' in ff[ikey]  :
				setattr(self,ikey,pd.read_hdf(filename,ikey)) 
				continue 
			# est-ce que c'est un group contenant uniquement d'autres groupes et/ou  :
			group_only = True 
			for isubkey in ff[ikey] : 
				ckey = ikey+'/'+isubkey
				if not type(ff[ckey]) == h5py._hl.group.Group  :
					group_only = False 
			# c'est un group contenant uniquement des autres group 'sta/FR' : on ne lit rien  
			if group_only == True : 
				continue 
			# ce group contient des datasets on lit sous forme de dict : 
			setattr(self,ikey,read_group_as_dict(ff[ikey]))
		# on a un dataset 
		else :
			setattr(self,ikey,ff[ikey][()])
	ff.close()

#=========================================================================================
#
# FUNCTION READING/WRITING TO AN H5 FROM ITS FILENAME, CALLING OTHER FUNCTION BELOW 
#
#============================================================================================
def filename_read_group_as_dict(filename,group_name) : 
	try :
		ff = h5py.File(filename,'r')
	except : 
		ff = h5py.File(filename,'a')

	return read_group_as_dict(ff[group_name])

def filename_write_dict_to_group(filename,dict_,prefix='/') : 
	try :
		ff = h5py.File(filename,'w')
	except : 
		ff = h5py.File(filename,'a')
		
	copy_dict_to_group(ff,dict_,prefix)
	ff.close()


#=========================================================================================
#
# FUNCTIONS WORKING ON AN ALREADY OPENED H5 FILE : WE PASS THE GROUP TO IT : 
#
#========================================================================================
def read_group_as_dict(group) :
	db = {}
	for ikey in group.keys() :	
		try :
			if group[ikey].ndim > 0 :
				db[ikey] = group[ikey][:]
			else :
				db[ikey] = group[ikey][()]
		except : 
			pass 
	return db


def copy_dataset_tree(h5_src,h5_dest,dname) : 
	''' copy : h5_src['/dname/sta1/sta2/*'] to h5_dest if to does not exist'''
	for sta1 in h5_src[dname] :
		for sta2  in h5_src[dname][sta1]:
			dset_root = dname+'/'+sta1+'/'+sta2
			if not dset_root in h5_dest : 
				h5_dest.copy(h5_src[dset_root],dset_root)


def copy_dict_to_group(h5,dict_,prefix=''):
	''' copy the dict dict_ into h5[/prefix/key_list] '''
	for ikey in dict_ :
		#	if ikey=='md' : ipdb.set_trace()
		h5.create_group(prefix+'/'+ikey)
		for dname, dset in dict_[ikey].items():
			#if dname=='id' : ipdb.set_trace()
			if type(dset)==list and len(dset) > 0 and type(dset[0]) == dict :
				copy_dict_to_group(h5,{dname : dset[0]},prefix='/'+ikey)
			if type(dset) == dict : 
				copy_dict_to_group(h5,{dname : dset},prefix='/'+ikey)

			# list of string ['ttoto','tata'] 
			if type(dset)==list and len(dset) > 0 and isinstance(dset[0],str):
				dset_tmp = []
				for mot in dset:
					dset_tmp.append(mot.encode())
				h5.create_dataset(prefix+'/'+ikey+'/'+dname,data=dset_tmp)
			# for matlab : replace False by 0 and True by 1 :
			else :
				if type(dset) == bool :
					if dset == False :
						dset = 0 
					if dset == True :
						dset = 1
				try :
					h5.create_dataset(prefix+'/'+ikey+'/'+dname,data=dset)
				except :
					pass 
	return h5 


				
