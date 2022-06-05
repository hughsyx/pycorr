################################################
# ava.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################

import pycorr.mods.s2d as s2d 
import pycorr.mods.cc2da as cc2da 
import pycorr.analyse.main_loop as main_loop 
import pycorr.mods.dd as dd 
import pycorr.mods.lang as lang 
import pycorr.analyse.main_loop as main_loop 
import numpy as np 
import time

from scipy.interpolate import interp1d 

try : 
	import ipdb
except :
	pass


def ava_ml(in_dir,in_ava,in_ml) : 
	in_={}
	in_['p1']=5              # used to filter + defining coda windows :
  	in_['p2']=10              	
	in_ava = lang.parse_options(in_,in_ava)
	tag='pydb_ava'+'__'+str(int(in_['p1'])).zfill(2)+'_'+str(int(in_['p2'])).zfill(2)+'s'
	in_ava['tag'] = tag 
	in_ml['load_cc'] = False 
	main_loop(in_dir,ava,in_ava,in_ml)

def ava(cc,in_ava) :
	# return null measurements if the surface waves could not be found on the CC
	in_={}
	in_['p1']=5              # used to filter + defining coda windows :
  	in_['p2']=10              	
	in_ = lang.parse_options(in_,in_ava)

	#preparing output : 
	out={}
	out['ava']=np.zeros(2)
	#.. 
	md_c={}
	md_c['date1']=cc['date1'][0]
	md_c['date2']=cc['date2'][-1]
	md_c['note']='ampl. computed on the ref correlation causal then acausal'

	ndate=len(cc['date1'])
	nwf  =len(cc['t'])
	dt = cc['t'][1]-cc['t'][0]
	ref=s2d.filter(cc['ref'],in_['p1'],in_['p2'],dt);
	win=cc2da.get_window(cc['ref'],cc['t'],cc['dist'],{'plot':'False','p1':in_['p1'],'p2':in_['p2']})
	if not win['success'] : 
		return out, md_c, in_ 

	I1=int(win['pos']['sw'][0])
	I2=int(win['pos']['sw'][1])
	#..
	I3=int(win['neg']['sw'][0])
	I4=int(win['neg']['sw'][1])
	#..	
	out['ava'][0] =	max(abs(cc['ref'][I1:I2]))
	out['ava'][1] =	max(abs(cc['ref'][I3:I4]))

	#cc has the following field : 
	# cc['date1'] : 
	# cc['date1'] :  numpy.ndarray  [ 735600.04166667  735600.08333333  735600.125  
	# cc['date2'] :  numpy.ndarray  [ 735600.08333333  735600.125       735600.16666
	# cc['dist']  :  float  0.0
	# cc['cc_mat']:  h5py._hl.dataset.Dataset  <HDF5 dataset "ZZ": shape (176, 6001), type "<f4
	# cc['t']     :  numpy.ndarray  [-600.  -599.8 -599.6 ...,  599.6  599.8  600. ]
	# cc['ref']   :  numpy.ndarray  [ 0.  0.  0. ...,  0.  0.  0.]

	# md_c : contains all fields common to all processed path such as dates  :
	# md_c={} 
	# md_c['date1']= cc['date1'][I1]
	# md_c['date2']= cc['date2'][I2]
	# md_c['date'] = np.mean((md_c['date1'],md_c['date2']),axis=0) #to be 

	# out : result of our measurements for a given path : 
	# out['dvv']=np.zeros((nstack))
	# out['coeffmatrix'] =np.zeros((nstack,nstretch))
	# out['corrcoef']   =np.zeros((nstack))
	#
	return out, md_c, in_ 