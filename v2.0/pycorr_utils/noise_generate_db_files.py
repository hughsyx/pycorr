################################################
# noise_generate_db_file.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import m_pycorr.mods.dd as dd
import glob 
import h5py
import os
import ipdb
import pickle
import sys 


def generate_db_files(in_dir='data_5.0hz/daily/FR/2019',station_list = './station.lst') :
	''' generate a db.pkl file into a data_1.0hz/daily/FR/2019 directory. This db is useful to compute CC with pycorr.
		This function is useful when using data not downloaded with pycorr'''

	dd.dispc('generating a db file for '+in_dir+' using'+station_list,'y','b')
	dd.dispc('  we assume that each day start at 00h00:00 i.e cut_len=9','g','n')
	dd.dispc('  and that data are stored in a directory whose name having the for data_1.0hz','g','n')
	dd.dispc('  from which we can extract the sampling rate','g','n')
	db = dict() 
	db['in_']= dict()
	db['in_']['tag'] = in_dir.split('/')[-2]
	db['in_']['pp']={}
	db['in_']['pp']['freq'] = float(in_dir.split('/')[0].split('_')[1][0:-2])
	db['in_']['pp']['cut_len']=0 

	db['sta']= dict()


	# if db_file already exist we rename it in order not to overwrite it :
	if os.path.isfile(in_dir+'/db.pkl') : 
		try :
			os.rename(in_dir+'/db.pkl',in_dir,'/db_bkp.pkl')
		except	:	
			pass 


	# read station list from file station.lst :
	sta_list = read_station_list(station_list)

	# loop on daily noise file : data_1.0hz/daily/FR/day*.h5 
	day_list = glob.glob(in_dir+'/*.h5') 
	for iday in day_list :
		#print(iday)
		# list station id in this day of data 
		id_list = h5_get_id_list(iday)
		for iid2 in id_list :
			iid = iid2.replace('.','_') # dp.pkl on a dies ID type : FR_ISO_00 et non FR_ISO.00 comme dans les day*h5.
			# we do not have this station in our db :
			if not iid in db['sta'] :
				db['sta'][iid]=dict()
				#print(iid)
				# do we have this station id in the station file :
				if iid in sta_list :
					# create new entry 
					db['sta'][iid]={}
					for ifield in sta_list[iid] :
						db['sta'][iid][ifield] = str(sta_list[iid][ifield])

	dd.dispc('  found '+str(len(db['sta']))+' stations','c','n')
	# write the new db_file
	filename  = in_dir+'/db.pkl' 
	ff        = open(filename,'wb')
	pickle.dump(db,ff)
	ff.close()


def read_station_list(filename) :
	sta={}
	ff=open(filename,'r')
	for iline in ff : 
		(dc,net,name,loc,lat,lon,elev,depth) = iline.split()
		if loc=='--' : loc = ''
		if loc == u'' : kname=net+'_'+name+'_00'
		else: kname=net+'_'+name+'_'+loc 
		sta[kname]={}
		sta[kname]['net']  = net
		sta[kname]['name'] = name
		sta[kname]['loc']  = loc
		sta[kname]['kname']= kname
		sta[kname]['lat']  = float(lat)
		sta[kname]['lon']  = float(lon)
		sta[kname]['elev'] = float(elev)
		sta[kname]['depth']= float(depth)
		sta[kname]['dc']   = dc
	ff.close()
	return sta 





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
			id_ = inet+'_'+sta
			id_list.append(id_)

	return id_list


	

if __name__=="__main__" :
		#if executed as a script : 
		#!/usr/bin/env python3 
		if len(sys.argv)>1 :
			generate_db_files(in_dir=sys.argv[1],station_list = sys.argv[2])
		else :
			dd.dispc('USAGE : noise_generate_db_files.py data_1.0hz/daily/FR/2019 ./station.lst','c','n')


