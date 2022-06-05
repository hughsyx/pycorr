################################################
# pcnoise_plot_station_map_single_day.py
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3

import sys, pickle, h5py, glob
import obspy.core as obs
import pycorr_mods.dd as dd 
import pycorr_mods.map as map 
import numpy as np 
import ipdb


def print_welcome_message()  : 
	print('')
	dd.dispc('pcnoise_plot_station_map_single_day.py','y','b')
	dd.dispc('  plot stations of a single day or a single year : ','y','n')
	dd.dispc('  - in red  : stations that are in the db.pkl file','y','n')
	dd.dispc('  - in green: stations that are in the day_xxx.h5 files','y','n')
	dd.dispc('  if data from all stations have been downloaded everything should be in green','y','n')

	print('')
	dd.dispc('  USAGE :','w','b')
	dd.dispc('  from a ipython session:','w','b')
	dd.dispc('  import pycorr_utils.pcnoise_plot_station_map_single_day as toto','w','n')
	dd.dispc('  toto.main(\'data_5.0hz/daily/CH/2016/\') ','w','n')
	dd.dispc('  toto.main(\'data_5.0hz/daily/CH/2016/day_002.h5\') ','w','n')
	print('')
	dd.dispc('  from the terminal:','w','b')
	dd.dispc('  pcnoise_plot_station_map_single_day.py data_5.0hz/daily/CH/2016/day_002.h5','w','n')
	dd.dispc('  pcnoise_plot_station_map_single_day.py data_5.0hz/daily/CH/2016/','w','n')


def main(path_) : 
	# clean path_ and look what can of path we did specify  :
	if path_[-1] == '/' : 
		path_ = path_[0:-1]
	lvl = len(path_.split('/'))

	# we ask for a specific day i.e path = data_5.0hz/daily/CH/2016/day_002.h5
	if lvl == 5 : 
		plot_station_map_single_day(path_)
	# we ask for a specific year i.e path = data_5.0hz/daily/CH/2016/
	if lvl ==4 : 
		plot_station_map_one_year(path_)


def plot_station_map_one_year(path_) : 
	db = plot_db_file(path_) 

	#open the hdf5 file an plot the corresponding station in green:
	day_list = glob.glob(path_+'/day_???.h5') 
	day_list.sort()

	db_plot = {} # station to be plotted 
	for iday in day_list:
		dd.dispc('  adding stations of '+iday,'c','n')
		db_plot =update_db_plot_with_station_of_this_day(iday,db,db_plot)
	map.stations_from_dict_of_stations(db_plot,c='g',s=100,marker='^',alpha=1,edgecolors='k',label='',text={}) 



def plot_station_map_single_day(path_) :
	print('')
	dd.dispc('  working on '+path_,'c','n')

	# open db file an plot it : 
	db = plot_db_file(path_) 

	#open the hdf5 file an plot the corresponding station in green:
	db_plot = {} 
	db_plot =update_db_plot_with_station_of_this_day(path_,db,db_plot)
	map.stations_from_dict_of_stations(db_plot,c='g',s=100,marker='^',alpha=1,edgecolors='k',label='',text={}) 


def update_db_plot_with_station_of_this_day(path_this_day,db,db_plot) :
	h5 = h5py.File(path_this_day,'r')
	for inet in h5 :
		if inet == '_metadata' : continue 
		for ista in h5[inet] : 
			ista_key = inet+'_'+ista
			ista_key = ista_key.replace('.','_')
			if ista_key in db['sta'] :
				db_plot[ista_key] = db['sta'][ista_key]
			else : 
				dd.dispc('  '+ista_key+' is in the .h5 file but not in the db.pkl file','r','r')
	h5.close()
	return db_plot 


def plot_db_file(path_) :
	db_path = '/'.join(path_.split('/')[0:4])+'/db.pkl'
	with open(db_path, 'rb') as f:
		db = pickle.load(f)
	map.plot_station_map_from_dict_of_stations(db_path,color='r',topo=False)

	return db


if __name__=="__main__" :
	narg = len(sys.argv)
	if narg == 2 :
		main(sys.argv[1])
	else : 
		print_welcome_message()