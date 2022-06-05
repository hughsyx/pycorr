from obspy.core import UTCDateTime
import scipy.io as io
import pickle
import m_pycorr.mods.dd as dd



def main(path_to_db,in_) :
	# import pycorr_scripts.noise01_edit_station_db as toto
	# 
	#in_ is a dict containing input parameters that will overriden in the db.pkl file ['in_'] dictionnary 
	# EXAMPLE 
	# in_ ={}
	# in_['day2'] = UTCDateTime(2017,1,1)
	# toto.main('data_5.0hz/daily/IV/2016/db.pkl',in_)                                                                 

	print(path_to_db)

	db=pickle.load( open( path_to_db, "rb" ) )

	dd.dispc('original settings found in '+ path_to_db,'y','b')
	dd.dd(db['in_'])

	dd.dispc('','y','n')
	dd.dispc('going to be replaced by ','y','b')
	dd.dd(in_)

	for ikey in in_ : 
		db['in_'][ikey]=in_[ikey]

	pickle.dump( db, open( path_to_db, "wb" ))


