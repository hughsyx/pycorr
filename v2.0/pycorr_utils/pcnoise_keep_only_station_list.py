################################################
# pcnoise_keep_only_station_list.py
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import sys, os, glob, h5py
import pycorr_mods.dd as dd 
import pycorr_mods.lang as lang 
import ipdb 
import pickle


def print_welcome_message() : 
	print('')
	dd.dispc('pcnoise keep only station','y','b') 
	dd.dispc(' - scan a directory [data_5.0hz/daily/IV/2016] containing daily noise record ','w','n')
	dd.dispc(' - generate a new directory [data_5.0hz_tag/[...]containing only a subset of the station','w','n')
	print('')
	dd.dispc('  USAGE','w','b')
	dd.dispc('  from the terminal : NONE ','w','b')
	print('')
	dd.dispc('  from a script:','w','b')
	dd.dispc('  import pycorr_utils.pcnoise_keep_only_station_list as pcn','w','n')
	dd.dispc('  db="data_5.0hz/daily/IV_nrca/2016"','w','d')
	dd.dispc("  sta_list =['/IV/NRCA.00','/IV/PESA.00']",'w','d')
	dd.dispc("  in_ = {}",'w','d')
	dd.dispc("  in_['tag'] ='around_nrca'",'w','d')
	dd.dispc("  pcn.main(db,sta_list,in_)",'w','n')	
	print('')


def main(path_,station_list,inu) :
	# for a given set, year
	# determine output dir 
	# determine output file 
	# copy only a few dset. 

	in_ = {}
	in_['tag'] ='station_selected'
	in_ = lang.parse_options(in_,inu)

	# create ouput directory :	
	aa = path_.split('/')
	out_dir = aa[0]+'_'+in_['tag']+'/'+aa[1]+'/'+aa[2]+'/'+aa[3]

	try :
		os.makedirs(out_dir)
	except	:
		pass 

	# loop on day of the original database :
	day_list = glob.glob(path_+'/day_???.h5')
	for iday in day_list :
		out_file = out_dir+'/'+iday.split('/')[-1]
		if os.path.isfile(out_file) :
			dd.dispc('  '+out_file+' already exist','r','r')
			continue 
		dd.dispc('  '+out_file,'c','n')
		h5_in  = h5py.File(iday,'r')
		h5_out = h5py.File(out_file,'w')

		# copy _metadata
		if '_metadata' in h5_in : 
			h5_in.copy('/_metadata',h5_out['/'])

		#copy each dataset 
		for ista in station_list : 
			cnet = '/'+ista.split('/')[1]
			if ista in h5_in :
				if cnet not in h5_out : 
					h5_out.create_group(cnet) 
				h5_in.copy(ista,h5_out[cnet]) 
		h5_in.close()
		h5_out.close()

	#regenerate the db_file	for this set/year (ie IV/2016)
	db_out = out_dir+'/db.pkl' 
	if os.path.isfile(db_out) == False :
		db_in = path_+'/db.pkl' 
		_generate_db_file(db_in,db_out,station_list)


def _generate_db_file(db_in,db_out,station_list) :
	if os.path.isfile(db_in) == False : 
		dd.dispc('  could not generate the db.pkl file','r','r')
		return
	fin =open(db_in,'rb')
	db_pkl = pickle.load(fin)
	db_pkl_new={}
	if 'in_' in db_pkl :
		db_pkl_new['in_'] = db_pkl['in_']
	if 'sta' in db_pkl : 
		db_pkl_new['sta'] = {}
		for ista in station_list : 
			csta = ista.replace('/','_').replace('.','_')
			csta = csta[1:]
			if csta in db_pkl['sta'] :
				db_pkl_new['sta'][csta]  = db_pkl['sta'][csta]
	fin.close()
	fout = open(db_out,'wb')
	pickle.dump(db_pkl_new, fout)
	fout.close()


if __name__=="__main__" :
	narg = len(sys.argv)
	if narg == 2 :
		main(sys.argv[1])
	else : 
		print_welcome_message()
