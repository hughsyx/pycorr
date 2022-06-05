################################################
# pcnoise_prepare_for_redownloading.py
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import h5py 
import os, glob, sys, shutil, pickle 
import pycorr_mods.dd as dd 
import ipdb
import numpy as np 
import pickle 

def check_input_argument(path_,min_percent) :
	display_welcome_message()

	if path_[-1] == '/' : 
		path_ = path_[0:-1]

	lvl = len(path_.split('/'))
	# we ask for a specific year i.e path = data_5.0hz/daily/CH/2016
	if path_.split('/')[-1].isdigit() : 
		generate_done_files(path_, min_percent) 
	# data_5hz 
	elif lvl  == 1 :
		net_list = glob.glob(path_+'/daily/*/')
		for inet in net_list : 
			year_list = glob.glob(inet +'/*/')
			for iyear in year_list : 
				generate_done_files(iyear, min_percent)
	# data_5hz/daily 
	elif lvl == 2 :
		net_list = glob.glob(path_+'/*/')
		for inet in net_list : 
			year_list = glob.glob(inet +'/*/')
			for iyear in year_list : 
				generate_done_files(iyear, min_percent)
	# data_5hz/daily/CH 
	elif lvl == 3 :
		year_list = glob.glob(path_ +'/*/')
		for iyear in year_list : 
			generate_done_files(iyear, min_percent)
	print('')


def display_welcome_message() : 
	print()
	dd.dispc('prepare for redownloading','y','b')
	dd.dispc('- delete all .cplt file','y','n')
	dd.dispc('- for each set/year compare the size of each file to largest file','y','n')
	dd.dispc('- for files larger than X\% of the largest one, we create a _day_XXX.done file so that they','y','n')
	dd.dispc('  won''t be redownloaded','y','n')
	dd.dispc('- for files smaller than X\% of the largest one, we make sure there is no _day_XXX.done file','y','n')
	dd.dispc('  so that we can redownload them','y','n')
	dd.dispc(' we can then run get_data in mode 1 to complete the dataset','y','n')

def	generate_done_files(path_, min_percent=90) :
	''' path : string that is the path to a specific net/year e. g : data_5.0hz/daily/CH/2016
		min_percent : [int in percent] we will re-download i.e not generate a done file, for all day of data if more than min_percent is missing. 
	'''
	print('')
	dd.dispc('  working on '+path_,'c','b')

	h5_day  = glob.glob(path_+'/day*h5')
	h5_size = np.zeros((len(h5_day)))

	for inum,iday in enumerate(h5_day) : 
		h5_size[inum] = os.path.getsize(iday)
	#print(h5_size)
	min_size = min_percent*h5_size.max()/100
	nday_large = len(h5_size[h5_size > min_size])
	nday_small = len(h5_size[h5_size <=  min_size]) 
	dd.dispc('  all files that are less than '+str(min_size/1000000)+' Mo will be downloaded','c','n')
	dd.dispc('  i.e we will take care that they do not have any .done file','c','n')
	dd.dispc('    '+str(nday_large)+' daily files are large enough and will have a done file','c','n')
	dd.dispc('    '+str(nday_small)+' daily files are too small and will not have a done file','c','n')

	for iday in h5_day : 
		done_file = path_+'/_'+iday.split('/')[-1][0:-3]+'.done' 
		cplt_file = iday[0:-3]+'.cplt'
		# fichier is large enough. It is probably complete => create a done file so that it is not downloaded
		if os.path.getsize(iday) > min_size :
			if not os.path.isfile(done_file) :
				create_lock_file(done_file)
		# file is too small : we should complete it => delete any previous done file
		else :
			if os.path.isfile(done_file):
				os.remove(done_file)
		# in any case we remove .cplt file : 
		if os.path.isfile(cplt_file) :
			dd.dispc('      removing '+cplt_file,'r','n')
			os.remove(cplt_file) 
	#	print(iday)

def create_lock_file(filename) : 
    ff=open(filename,'wb')
    pickle.dump(filename,ff)
    ff.close()



if __name__=="__main__" :
	#if executed as a script : 
	if len(sys.argv)==2 :
		check_input_argument(sys.argv[1],min_percent=90) 
		#generate_done_files(sys.argv[1],90)
	elif len(sys.argv)==3 : 
		check_input_argument(sys.argv[1],min_percent=float(sys.argv[2]))
		#generate_done_files(sys.argv[1],min_percent=float(sys.argv[2]))
	else :
		display_welcome_message()
		print()
		dd.dispc('USAGE : noise_generate_done_files data_5.0hz/daily/CH/2016','w','b') 
		dd.dispc('      : noise_generate_done_files data_5.0hz/daily/CH','w','b') 
		dd.dispc('      : noise_generate_done_files data_5.0hz/daily','w','b') 
		dd.dispc('      : noise_generate_done_files data_5.0hz','w','b')  
		dd.dispc('      : noise_generate_done_files data_5.0hz/daily/CH/2016 90','w','b')
		dd.dispc('           generate done file for file that are larger than 90\% of the largest file','w','n') 
		dd.dispc('      : noise_generate_done_files data_5.0hz/daily/CH/2016 100','w','b') 
		dd.dispc('           re-dl all files smaller than the largest one','w','n')
		dd.dispc('           =>remove all .done and .cplt files','w','n') 
		dd.dispc('      : noise_generate_done_files data_5.0hz/daily/CH/2016 0','w','b') 
		dd.dispc('           create .done files for all daily file, and remove all .cplt file','w','n') 
		print()