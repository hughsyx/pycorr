import numpy as np, pandas as pd 
import h5py , ipdb 

from m_pycorr_dev.mods import dd 


# TODO : ajouter une liste de fichier dans db_c1_sta.h5 
# depend de db_c1_sta, db_c1*


# convert input to indice : voir commment packer/unpacker des variables : 
# def _read_c1_ref()
# def _ref_c1_matrix()
#
# on fait une classe uniquement pour ne pas faire des acces disques a chaque fois qu'on cherche des infos 
# sur les stations => suppose que la liste des stations tient en RAM

class C1_read : 
	def __init__(self,C1_dir) : 
		self.c1_dir = C1_dir 
		self.db    = self.c1_dir+'/db_c1.h5' 
		self.sta   = pd.read_hdf(self.db,'sta')  # tableau unique toute les stations
		#self.file  = pd.read_hdf(self.db,'file') # useless 
		self.cmp   = pd.read_hdf(self.db,'cmp')  # useless


	def from_file(self,kfile,I_path,I_cmp='ZZ',type_='ref',t1=0,t2=np.inf) :
		# convert kfile from str to int if needed 
		# convert argument from str/scalar to list if needed : => refactored later :
		#if type(I_path) == str : print(I_path)    # to be done later 
		#if type(I_path) == int : I_path =[I_path]
		db = self.c1_dir+'/db_c1_by_file.h5'

		I_path, I_cmp = self._read_from_parse_input(I_path,I_cmp)

		# get the filename corresponding to kfile (if the user did not put directly the filename into kfile) 
		# normally we should have used the db_c1_from_files_tree 
		if not type(kfile) == str :
			file_table = pd.read_hdf(self.db,'file',where=[kfile]) # contain both filename + ncpl 
			filename = file_table.index[0]
		else : 
			filename  = kfile 
		filename_full = self.c1_dir+'/'+filename+'.h5'

		# initialize list of station indice for each path requested by user : 
		cpl = self._cpl_init(I_path,db,filename) 
		self._cpl_complete_metadata_from_db_c1(cpl)
		return self._read_c1_ref_from_single_file(filename_full,cpl,I_cmp) 


	def _read_c1_ref_from_single_file(self,filename,cpl,kcmp) : 
		''' low level function that reads several reference C1 inside the same file, 
			for a given compononent and list of couple 

			INPUTS 
			cpl : dictionnary listing the couple of station to be read  : 
				cpl['I_sta'][kpath,:]   : 2D matrix of couple of station inside a given file to be read
				cpl['dist']  = 1D vector of inter-station distance 
				cpl['az']    = 1D vector of azimuth 
				cpl['baz']   = 1D vector of baz 
				cpl['filename'] : name of the file from which we read the correlation 
			kcmp : string of the component name to be read : "ZZ","ZN"
			'''
		cc = {} # dictionary containing the output : 
		cc['cpl'] = cpl 

		ff = h5py.File(filename,'r')
		ncpl = len(cpl['I_sta']) 
		for icpl in range(0,ncpl) : 
			sta0 = self.sta.index[cpl['I_sta'][icpl][0]]
			sta1 = self.sta.index[cpl['I_sta'][icpl][1]]
			dset_name = '/ref/'+sta0+'/' + sta1+'/'+kcmp
			if icpl == 0 : 
				trace = ff[dset_name][()]
				nwf = len(trace)
				cc['tr'] = np.zeros((ncpl,nwf))
			cc['tr'][icpl,:] = ff[dset_name][()]
		ff.close()
		return cc 


	def _cpl_init(self,I_path,db,dset_name) :
		''' initialize list of station indice for each path requested by user from a db_c1_from_files.h5 file'''
		ff = h5py.File(db,'r')
		cpl = {}
		cpl['I_sta']= ff['/'+dset_name+'/I_sta'][I_path,:]
		cpl['dist'] = ff['/'+dset_name+'/dist'][I_path]
		cpl['az']   = ff['/'+dset_name+'/az'][I_path]
		cpl['baz']  = ff['/'+dset_name+'/baz'][I_path]
		ff.close()
		return cpl 

	def _cpl_complete_metadata_from_db_c1(self,cpl) :
		''' complete cpl dictionnaries, adding lon,lat,elev,depth from db_c1.h5 station table using cpl[I_sta] = 2D matrix of station indice'''
		cpl['lon']   = self.sta['lon'][cpl['I_sta']]
		cpl['lat']   = self.sta['lat'][cpl['I_sta']]
		cpl['elev']  = self.sta['elev'][cpl['I_sta']]
		cpl['depth'] = self.sta['depth'][cpl['I_sta']]


	def _read_from_parse_input(self,I_path,I_cmp) :
		''' parse input of read_from* function. Convert I_path to list, and I_cmp to str : '''
		# convert I_path into a list if it is a int : 
		if type(I_path) == int : 
			I_path= [I_path]
		# convert I_cmp to a string "ZZ" if it is an component indice 
		if type(I_cmp) == int : 
			I_cmp = self.cmp.iat[0,I_cmp] # MAYBE WRONG
		return I_path, I_cmp



#==================================================================




#def _parse_input() :
#	# kfile string => indice 
#	# I_path = 'FR.OGAG.00/FR.RSF.00'-> nothing
#	# I_cmp : 'ZZ' -> 1 
#	pass 


#def _get_filename_from_indice(C1_dir,I_file) : 
#	return pd.read_hdf(C1_dir+'/db_c1.h5','file',where = [I_file]).index[0]+'.h5'
	
#def _get_station_pair_metadata_from_indice(C1_dir,I_sta) :
#	return  pd.read_hdf(C1_dir +'/db_c1.h5' ,'sta',where=I_sta)



	# fonction de bas niveau :
	# 2. initialiser matrice de CC 
	# 3. boucle sur les station 
	#   -- construit le dset_name 
	#   -- lit les donnees
	#   -- renvoit le resultat 
