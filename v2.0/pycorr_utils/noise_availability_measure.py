################################################
# noise_availability.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import glob , pickle
import numpy as np
import h5py, pickle
import os, sys
from obspy.core.utcdatetime import UTCDateTime
from scipy.io import savemat 


#try : 
#	import ipdb
#except 	:	
#	pass 


def measure_station_availability_for_whole_db(in_dir='./data_1.0hz') :
	''' measure station availability for a whole db 
		OUTPUT : a matlab file in availability/av_setname.mat indicating for each year/day/station/
		if data are available.

		Can be run in // : each process will then work on a separate set 

		RESULTS can be then plotted using pycorr.utils.noise_availability_plot

	'''

	net_list = glob.glob(in_dir+'/daily/*/');
	for inet in net_list :
		measure_station_availability_for_one_set(inet) 


def measure_station_availability_for_one_set(in_dir='./data_1.0hz/daily/TEST/'):
	''' measure station availability for a single set all year
		OUTPUT : a matlab file in availability/av_setname.mat indicating for each year/day/station/
		if data are available 

		RESULTS can be then plotted using pycorr.utils.noise_availability_plot
	'''

	# define output file name : 
	out_dir = 'availability' 
	if not os.path.isdir(out_dir) :
		os.mkdir(out_dir) 

	if in_dir[-1]=='/' : 
		in_dir = in_dir[:-1]

	# determining output name :
	out_root   = out_dir +'/av_'+ in_dir.split('/')[-1]
	out        = out_root+'.mat'
	out_lock   = out_root+'.'+str(os.getpid())+'.lock' 
	out_search = out_root+'.*' 		

	# check if it is already (being) done : 
	if len(glob.glob(out_search)) > 0 :
		dispc('    '+out_search+'already exist => returning','r','b')
		return
		
	create_lock_file(out_lock)
	dispc('    result will be written in '+ out,'c','b')

	# av will contain all information about station avail :
	av={}
	av['sta'] = {}

	# first list all year : 
	year_dir_list=glob.glob(in_dir+'/*/')
	year_dir_list.sort()
	year1 = int(year_dir_list[0].split('/')[-2])
	year2 = int(year_dir_list[-1].split('/')[-2])

	# now construct the matlab compatible datevector 
	day1 = UTCDateTime(year1,1,1).toordinal()
	day2 = UTCDateTime(year2,12,31).toordinal()
	av['date'] = np.arange(day1,day2+1)
	av['date'] = av['date']+366 # for matlab compatible

	# construct the python compatible vector : 
	day1 = UTCDateTime(year1,1,1).timestamp
	day2 = UTCDateTime(year2,12,31).timestamp
	av['datev'] = np.arange(day1,day2+1,86400)
	av['datev'] = av['datev'] 
	ndate = len(av['date'])

	# read station list from db_file : 
	db =read_and_merge_db(in_dir,year1,year2)

	# now loop on each day 
	kday=-1 # day index 
	year0 = 0 
	for iday in av['datev'] :
		kday = kday + 1
		cdate = UTCDateTime(iday)
		if cdate.year != year0  :
			year0 = cdate.year 
			dispc('    '+str(year0),'c','n')
		year  = str(cdate.year) 
		cday  = cdate.julday
		path = in_dir+'/'+year+'/day_'+str(cday).zfill(3)+'.h5'
		if not os.path.isfile(path) : 
			continue 
		id_list =  h5_get_id_list(path)
		for iid in id_list :
			if iid not in av['sta'].keys() :
				av['sta'][iid]=np.zeros(ndate)
			av['sta'][iid][kday]= 1 

	# convert av into a matlab friendly file
	mat={}
	mat['in_dir']= in_dir 
	mat['date']=av['date']
	mat['date_python']=av['datev']
	mat['id'] = [*av['sta'].keys()]

	nsta = len(mat['id'])
	nday = len(mat['date'])
	
	mat['lon'] = np.zeros((nsta))
	mat['lat'] = np.zeros((nsta))

	mat['av']=np.zeros((nsta,nday)); 
	ksta=-1 
	for ista in mat['id'] :
		ksta=ksta+1
		mat['av'][ksta,:] = av['sta'][ista]
		mat['lon'][ksta]  = db['sta'][ista.replace('.','_')]['lon']
		mat['lat'][ksta]  = db['sta'][ista.replace('.','_')]['lat']

	savemat(out,mat)

	#delete tmp file:
	delete_lock_file(out_lock)

def	read_and_merge_db(in_dir,year1,year2) :
	db_dest={}
	db_dest['sta']={}

	for iyear in range(year1,year2+1) :
		path_db_source = in_dir+'/'+str(iyear)+'/db.pkl'	
		if not os.path.isfile(path_db_source) : 
			continue 
		with open(path_db_source, 'rb') as fsrc:
			db_src  = pickle.load(fsrc)
			db_dest['sta']={**db_src['sta'], **db_dest['sta']}
	return db_dest 


def h5_get_id_list(path) : 
	id_list = []
	try :
		ff = h5py.File(path,'r')
	except : 
		return []
	net_list = ff.keys()
	for inet in net_list :
		if inet == '_metadata' : 
			continue
		for sta in ff[inet].keys() :
			id_ = inet+'.'+sta
			id_list.append(id_)

	return id_list



def delete_lock_file(fname) :
	if os.path.isfile(fname) : 
		try : 
			os.remove(fname) 
		except :
			pass 



def create_lock_file(fname) : 
	a=0 ; 
	with open(fname,'wb') as f : 
		pickle.dump(a,f) 



def dispc(text,color,attr) :
    # print text in color with attribute attr ! 
    if color =='r' :
        ncolor=31
    elif color=='g' :
        ncolor=32 
    elif color=='y' :
        ncolor=33
    elif color=='b' :
        ncolor=34
    elif color=='m' :
        ncolor=35
    elif color=='c' :
        ncolor=36
    elif color=='w' :
        ncolor=37
    elif color=='gray' :
        ncolor=38

    if attr=='b' :
        prefix='\x1b[1m';
    elif attr=='d' :
        prefix='\x1b[2m';
    elif attr=='u' :
        prefix ='\x1b[4m';
    elif attr=='blink' :
        prefix ='\x1b[5m'
    elif attr=='r' :
        prefix='\x1b[7m'
    else :
        prefix=''

    str_=prefix+'\x1b['+str(ncolor)+'m'+text+'\x1b[0m';
    print(str_)



# if executed as a script :
if __name__=="__main__" :
		if len(sys.argv)>1 :
			in_dir = sys.argv[1]
			if len(in_dir.split('/')) > 1 :
				measure_station_availability_for_one_set(in_dir=in_dir)
			else :
				measure_station_availability_for_whole_db(in_dir=in_dir)
		else :
			dispc('USAGE : get_station_availability data_1.0hz/daily/FR/ ','c','n')
