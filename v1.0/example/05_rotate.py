import m_pycorr.m40_rotate as rotate


in_ = {}
in_['path']      = 'C1_04_xcorr_cifalps_conti__xcorr_norm/'
in_['mode']      = 'ref_only' # 'all' 
in_['cc_cmp']    = ['ZZ','RR','TT'] #...
in_['format']    = 'float16' #...


rotate.main_rot(in_)