################################################
# noise_rename_set_starting_with_number
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import sys, h5py
import glob
import ipdb
import os
import m_pycorr.mods.dd as dd

def rename_set_starting_with_number(in_dir='data_1.0hz') :
	set_list = glob.glob(in_dir+'/daily/*/') 
	for iset in set_list : 
		iset_name = iset.split('/')[-2] 
		if iset_name[0].isdigit() :
			iset_newname = 'x'+iset_name 
			dest =  '/'.join(iset.split('/')[0:-2])+'/'+iset_newname
			dd.dispc('  '+iset+' -> '+dest,'c','n')
			os.rename(iset,dest)

if __name__ =="__main__" :
	if len(sys.argv) == 2:
		print(sys.argv)
		rename_set_starting_with_number(in_dir=sys.arg[1]) 
	elif len(sys.argv) == 1 : 
		datadir_list = glob.glob('data_*hz') 
		for idatadir in datadir_list :
			rename_set_starting_with_number(in_dir = idatadir)
	else  :
		print('')
		print('USAGE : noise_rename_set_starting_with_number data_1.0hz/daily')
		print('   OR : noise_rename_set_starting_with_number')
		print('')
	#noise.print_number_of_days(net='CH',year=2016,in_dir='data_5.0hz') 
