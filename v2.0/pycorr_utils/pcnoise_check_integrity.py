################################################
# pcnoise_check_integrity
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import sys, os, glob, h5py
from datetime import datetime
import pycorr_mods.dd as dd 

def print_welcome_message() : 
	print('')
	dd.dispc('pcnoise_check_integrity.py','y','b')
	dd.dispc('check the integrity of daily noise records files day_xxx.h5','y','n')
	dd.dispc('i.e make sure that all trace can be read, and that there are _metadata in each file','y','n')
	dd.dispc('It creates a _check_integrity directory that contains two list of stations :','y','n')
	dd.dispc('for instance if we have a set CH and for the year 2016 : ','y','n')
	dd.dispc('  - CH_2016_fine.txt : list all files that are fine','y','n')
	dd.dispc('  - CH_2016_wrong.txt : list all files that are wrong','y','n')
	dd.dispc('later wrong file can be deleted and redownloaded using pcnoise_prepare_for_redownloading.py','y','n')
	print('')
	dd.dispc('USAGE','w','b')
	dd.dispc('from the terminal:','w','b')
	dd.dispc('check_integrity.py data_5.0hz/daily/CH/2016/','w','n')
	dd.dispc('check_integrity.py data_5.0hz/daily/CH/','w','n')
	dd.dispc('check_integrity.py data_5.0hz/daily/','w','n')
	dd.dispc('check_integrity.py data_5.0hz','w','n')
	print('')
	dd.dispc('from a script:','w','b')
	dd.dispc('import pycorr_utils.pcnoise_check_integrity as toto','w','n')
	dd.dispc('toto.main(\'data_5.0hz/daily/CH/2016/\')','w','n')
	print('')

def main(path_) :

	# create output directory 
	now  = datetime.now()
	date_= now.strftime('%Y_%m_%d__%Hh%M') #'S'
	
	out_dir = '_check_integrity_'+date_ 
	try :
		os.mkdir(out_dir)
	except : 
		pass 

	# clean path_ and look what can of path we did specify  :
	if path_[-1] == '/' : 
		path_ = path_[0:-1]
	lvl = len(path_.split('/'))

	# we ask for a specific year i.e path = data_5.0hz/daily/CH/2016
	if path_.split('/')[-1].isdigit() : 
		check_integrity(path_, out_dir) 
	# data_5hz 
	elif lvl  == 1 :
		net_list = glob.glob(path_+'/daily/*/')
		for inet in net_list : 
			year_list = glob.glob(inet +'/*/')
			for iyear in year_list : 
				check_integrity(iyear, out_dir) 
	# data_5hz/daily 
	elif lvl == 2 :
		net_list = glob.glob(path_+'/*/')
		for inet in net_list : 
			year_list = glob.glob(inet +'/*/')
			for iyear in year_list : 
				check_integrity(iyear, out_dir) 
	# data_5hz/daily/CH 
	elif lvl == 3 :
		year_list = glob.glob(path_ +'/*/')
		for iyear in year_list : 
			check_integrity(iyear, out_dir) 

def check_integrity(path_,out_dir) : 
	print('')
	if path_[-1] == '/' : 
		path_ = path_[0:-1]
	dd.dispc('  working on '+path_,'c','b')

	# create the output file that will contain the list of file which are not ok 
	year = path_.split('/')[-1] 
	set  = path_.split('/')[-2]
	out_file_wrong = out_dir+'/'+set+'_'+year+'_wrong.txt'
	out_file_fine  = out_dir+'/'+set+'_'+year+'_fine.txt'
	file_wrong = open(out_file_wrong,"w") 
	file_fine  = open(out_file_fine,"w") 

	h5_day  = glob.glob(path_+'/day*h5')
	h5_day.sort()

	for iday in h5_day : 
		is_ok, is_metadata = check_this_daily_noise_record(iday)
		if is_ok== True and is_metadata==True:
			dd.dispc('  '+iday+' integrity = '+str(is_ok)+'  is_metadata = '+str(is_metadata),'c','n')		
			file_fine.write(path_+'/'+iday+'\n')
		else :
			dd.dispc('  '+iday+' integrity = '+str(is_ok)+'  is_metadata = '+str(is_metadata),'r','r')		
			file_wrong.write(path_+'/'+iday+'\n')
	file_wrong.close()
	file_fine.close()


def check_this_daily_noise_record(kday) :
	is_metadata = False 
	try :
		h5_day = h5py.File(kday,'r') 
		for inet in h5_day : 
			if inet == '_metadata' : 
				is_metadata=True
			else : 
				for ista in h5_day[inet] : 
					for ich in h5_day[inet][ista] :
						trace = h5_day[inet][ista][ich]
		return True, is_metadata 
	except:
		return False, False 



if __name__=="__main__" :
	narg = len(sys.argv)
	if narg == 2 :
		main(sys.argv[1])
	else : 
		print_welcome_message()
