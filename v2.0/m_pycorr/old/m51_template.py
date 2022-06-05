
import pycorr_x.modules.s2d as s2d 
import pycorr_x.modules.cc2da as cc2da 
import pycorr_x.analyse.main_loop as main_loop 
import pycorr.modules.dd as dd 
import pycorr.modules.lang as lang 
import pycorr_x.analyse.main_loop as main_loop 
import numpy as np 
import time

from scipy.interpolate import interp1d 

try : 
	import ipdb
except :
	pass


def xx_ml(in_dir,in_xx,in_ml) : 
	in_={}
	in_['p1']=5              # used to filter + defining coda windows :
  	in_['p2']=10              	
	in_xx = lang.parse_options(in_,in_xx)
	tag='pydb_dt'+'__'+str(int(in_['p1'])).zfill(2)+'_'+str(int(in_['p2'])).zfill(2)+'s'
	in_['tag'] = tag 
	main_loop(in_dir,dt,in_xx,in_ml)

def xx(cc,in_xx) :
	in_={}
	in_['p1']=5              # used to filter + defining coda windows :
  	in_['p2']=10              	
	in_ = lang.parse_options(in_,in_xx)

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