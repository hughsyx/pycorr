################################################
# pcnoise_generate_dbmat_from_dbpkl.py
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import sys, os, glob, pickle ,ipdb
import scipy.io as io
from datetime import datetime
import pycorr_mods.dd as dd 


def print_welcome_message()  : 
	print('')
	dd.dispc('pcnoise_generate_dbmat_from_dbpkl.py','y','b')
	dd.dispc('  read a data_5.0hz/daily/toto/2016/db.pkl file : ','y','n')
	dd.dispc('  and convert it to a db.mat if it does not already exist','y','n')
	dd.dispc('  useful to fix a data directory once it has been manipulated (merging dataset,...)','y','n')

	print('')
	dd.dispc('  USAGE :','w','b')
	dd.dispc('  from a ipython session:','w','b')
	dd.dispc('  import pycorr_utils.pcnoise_generate_dbmat_from_dbpkl as pcn ','w','n')
	dd.dispc('  pcn.main(\'data_5.0hz/daily/toto/2016\') ','w','n')
	print('')
	dd.dispc('  from the terminal:','w','b')
	dd.dispc('  pcnoise_generate_dbmat_from_dbpkl.py data_5.0hz/daily/toto/2016','w','n')


def main(path_) :
	now  = datetime.now()
	date_= now.strftime('%Y_%m_%d__%Hh%M') #'S'
	dd.dispc('  working on '+path_,'c','b')

	db_in = path_+'/db.pkl'
	db_out= path_+'/db.mat'
	if not os.path.isfile(db_in) :
		dd.dispc(db_in+' not found => exiting','r','r')
		return
	if os.path.isfile(db_out) :	
		dd.dispc(db_out+' already exist','r','r')
		return 

	dd.dispc('   reading '+db_in,'c','n')
	ff=open(db_in,'rb');
	db_ori=pickle.load(ff)
	ff.close()
	dd.dispc('   found '+str(len(db_ori['sta']))+' stations','c','n')

    #output db_file in a matlab-friendly way (taken from m01_get_data): 
	db = dict() 
	if 'ev' in db_ori :
		db['ev'] = db_ori['ev']
	db['in_']              = db_ori['in_']
	db['in_']['day1']      = float(db_ori['in_']['day1'].toordinal()) + 366.
	db['in_']['day2']      = float(db_ori['in_']['day2'].toordinal()) + 366.
	if db_ori['in_']['password'] is None:db['in_']['password']='None'
	if db_ori['in_']['user_id'] is None:db['in_']['user_id']='None'
	if db_ori['in_']['ev']['time_after'] is None:db['in_']['ev']['time_after']='None'
	db['sta']              = dict()
	db['sta']['lon']       = []
	db['sta']['lat']       = []
	db['sta']['sta']       = []
	db['sta']['net']       = []
	db['sta']['dc']        = []
	db['sta']['loc']       = []
	db['sta']['elev']      = []
	db['sta']['depth']     = []

	for ista in db_ori['sta'] :
		db['sta']['sta'].append(ista) 
		db['sta']['net'].append(db_ori['sta'][ista]['net'])
		db['sta']['dc'].append(db_ori['sta'][ista]['dc'])
		db['sta']['lon'].append(float(db_ori['sta'][ista]['lon']))
		db['sta']['lat'].append(float(db_ori['sta'][ista]['lat']))
		db['sta']['loc'].append(db_ori['sta'][ista]['loc'])
		db['sta']['elev'].append(float(db_ori['sta'][ista]['elev']))
		db['sta']['depth'].append(float(db_ori['sta'][ista]['depth']))

	dd.dispc('   writing '+db_out,'c','n')
	io.savemat(db_out,db)

if __name__=="__main__" :
	narg = len(sys.argv)
	if narg == 2 :
		main(sys.argv[1])
	else : 
		print_welcome_message()
