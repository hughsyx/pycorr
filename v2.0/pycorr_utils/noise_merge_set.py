################################################
# noise_merge_set.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import h5py 
import os, glob, sys, shutil, pickle 
import m_pycorr.mods.dd as dd 
import ipdb

def merge_set(path_source, path_dest,year_list) : 
	''' add all station data from path_source into path_dest for year in year list
		EXAMPLE : 
			import utils.pycorr_bin.noise_merge_set as nms
			path_source=[] 
			path_source.append('data_5.0hz/daily/CH/')
			path_source.append('data_5.0hz/daily/IT/')
			path_dest = 'data_5.0hz/daily/toto/'
			nms.merge_set(path_source,path_dest,2016)
	'''

	# make sure input are list :
	if not (type(year_list) is list) : 
		year_list = [year_list]

	if not type(path_source) is list : 
		path_source = [path_source]


	dd.dispc('-----','y','b')
	for k, ipath in enumerate(path_source): 
		if k==0 :
			dd.dispc('merging '+ipath,'y','b') 
		else 	:
			dd.dispc('        '+ipath,'y','b') 
	dd.dispc('into '+path_dest,'y','b')
	dd.dispc('-----','y','b')

	# loop on year :
	for iyear in year_list : 
		dd.dispc(' working on year '+str(iyear),'c','b')

		# create the output directory if necessary: 
		out_dir = path_dest + '/'+str(iyear)
		if not os.path.isdir(out_dir):
			os.makedirs(out_dir)

		# merge the db_file :
		merge_db_files(path_source, path_dest, iyear)

		# loop on day and copy new station into the dest :
		for iday in range(1,367):
			dd.dispc(' working on '+str(iyear) + '  '+str(iday).zfill(3),'c','n')
			suffix ='/'+str(iyear)+'/day_'+str(iday).zfill(3)+'.h5'
			for ipath in path_source :
				add_source_into_dest(ipath+suffix, path_dest+suffix)
				add_source_into_dest(ipath+suffix, path_dest+suffix)



def merge_db_files(path_source , path_dest, iyear ) : 

	# merge the db of the other source to the dest :
	db_dest={}
	db_dest['sta']={}
	for isource in path_source[0:] :
		path_db_source = isource+'/'+str(iyear)+'/db.pkl'	
		if not os.path.isfile(path_db_source) : 
			continue 
		with open(path_db_source, 'rb') as fsrc:
			db_src  = pickle.load(fsrc)
			db_dest['sta']={**db_src['sta'], **db_dest['sta']}
	db_dest['in_']=db_src['in_']
	db_dest['ev'] = db_src['ev']

	# remove the db.pkl of the dest folder if it exists
	path_db_dest = path_dest+'/'+str(iyear)+'/db.pkl'; 
	if os.path.isfile(path_db_dest) : 
		os.remove(path_db_dest) 

	# write into the dest
	with open(path_dest+'/'+str(iyear)+'/db.pkl','ab') as fdest :
		pickle.dump(db_dest,fdest)



def add_source_into_dest(file_source,file_dest):
	if os.path.isfile(file_source) :
		h5_source = h5py.File(file_source,'r')
		h5_dest   = h5py.File(file_dest,'a')

		if not '_metadata' in h5_dest : 
			try :
				h5_source.copy('/_metadata',h5_dest['/'])
			except :
				pass

		for inet in h5_source :
			if inet == '_metadata' : 
				continue
			if inet not in h5_dest : 
				h5_dest.create_group('/'+inet) 
				for ista in h5_source['/'+inet] :
					group_name = '/'+inet+'/'+ista 
					print(group_name)
					if group_name not in h5_dest :
						h5_source.copy(group_name,h5_dest['/'+inet]) 

		h5_source.close()
		h5_dest.close()


if __name__=="__main__" :
	print('x')
	#if executed as a script : 
	#!/usr/bin/env python3 
	if len(sys.argv)>3 :
		merge_set(sys.argv[1],sys.argv[2],sys.argv[3]) 
	else :
		dd.dispc('USAGE : noise_merge_set.py data_5.0hz/daily/CH/ data_5.0hz/daily/toto/ 2016','w','n')

