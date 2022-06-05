################################################
# dvv.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


import m_pycorr.mods.s2d as s2d 
import m_pycorr.analyse.mods.cc2da as cc2da 
import pycorr_mods.dd as dd 
import pycorr_mods.lang as lang 
import m_pycorr.analyse.main_loop as main_loop 
import numpy as np 
import time
from scipy.interpolate import interp1d 
from scipy.io import savemat
from scipy.fftpack import fft, ifft

try : 
	import ipdb
except :
	pass

def dvv_ml(in_dir,in_dvv,in_ml) : 
	in_={}
	in_['p1']     = 1.
	in_['p2']     = 3.
	in_['v1']     = 1.5
	in_['v2']     = 5.
	in_['coda_dp']= 5.
	in_['tf']     = 120.  #length of coda 
	in_['eqspec'] = False 
	in_['eqspec_ac'] = False # equalize spectrum separately on the causal/acausal
	in_['eqspec_full'] = False
	in_['stack']  = 10
	in_['stack_shift']= 5.
	in_['stack_norm']=False 
	in_['stretch']=np.arange(-2,2,0.1)*1e-3
	in_['export_matlab']=False
	in_['keep_coeffmatrix']=False
	in_['print_dvv']=False

	in_dvv = lang.parse_options(in_,in_dvv)
	tag='pydb_dvv'+'__'+str(int(in_['p1'])).zfill(2)+'_'+str(int(in_['p2'])).zfill(2)+'s'
	tag=tag+'__stack_'+str(in_['stack'])+'__'+'stack_shift_'+str(in_['stack_shift'])
	tag =tag+'__v1_'+str(in_['v1']) +'_v2_'+str(in_['v2'])+'_coda_dp_'+str(in_['coda_dp'])
	tag=tag+'__tf_'+str(int(in_['tf']))
	if in_['eqspec'] :
		tag=tag+'__eqspec'
	if in_['eqspec_ac'] : 
		tag=tag+'__eqspecac'
	if in_['eqspec_full'] : 
		tag=tag+'__eqspecFull'

	in_dvv['tag']=tag;
	main_loop.main_loop(in_dir,dvv,in_dvv,in_ml)

def dvv(cc,varargin):
	in_={}
	in_['p1']     = 1.
	in_['p2']     = 3.
	in_['v1']     = 1.5
	in_['v2']     = 5.
	in_['coda_dp']= 5.
	in_['tf']     = 120.  #length of coda 
	in_['eqspec'] = False 
	in_['eqspec_ac'] = False
	in_['eqspec_full'] = False
	in_['stack']  = 10
	in_['stack_shift']= 5.
	in_['stack_norm'] = False # normalize matrix of CC before stacking it ?
	in_['stretch']=np.arange(-2,2,0.1)*1e-3
	in_['export_matlab']= False
	in_['keep_coeffmatrix']=False
	in_['print_dvv']=False
        
	in_=lang.parse_options(in_,varargin)
	ndate=len(cc['date1'])
	nwf  =len(cc['t'])
	nstretch=len(in_['stretch'])
	dt = cc['t'][1]-cc['t'][0]

    #dd.dispc(str(cc['dist']),'y','b')
	#filtering reference : 	
	ref=s2d.filter(cc['ref'],in_['p1'],in_['p2'],dt);
	#determining coda window : we keep coda_length seconds of coda 
	in_win = {'p1':in_['p1'],'p2':in_['p2'],'sw_v1':in_['v1'],'sw_v2':in_['v2'],'coda_dp':in_['coda_dp']}
	win=cc2da.get_window(cc['ref'],cc['t'],cc['dist'],in_win)
	I_time1 = max((  1,np.round(win['neg']['coda'][1]-in_['tf']/dt)))
	I_time2 = min((nwf,np.round(win['pos']['coda'][0]+in_['tf']/dt)))

	I_win = np.concatenate((np.arange(I_time1,win['neg']['coda'][1]),np.arange(win['pos']['coda'][0],I_time2)))
	I_win=I_win.astype(int)
	#stretching reference : use 100 points of margin for dummy reasons
	dt_ext  =int(10/dt)
	dt_ext2 =int(5/dt)
	I_win_ext_start=int(max(0  ,I_time1-dt_ext))
	I_win_ext_end  =int(min(nwf,I_time2+dt_ext))
	I_win_ext = np.concatenate((np.arange(I_win_ext_start,int(win['neg']['coda'][1])),np.arange(int(win['pos']['coda'][0]),I_win_ext_end)))
	refs=np.zeros((nstretch,nwf))
	cc['t']=np.squeeze(cc['t'])
	f_interp = interp1d(cc['t'][I_win_ext],ref[I_win_ext],axis=0,kind='linear')
	##
	k=0
	for istretch in in_['stretch'] :
		timei=cc['t']*(1+istretch)
		refs[k,I_win_ext[0]+dt_ext2:I_win_ext[-1]-dt_ext2]=f_interp(timei[I_win_ext[0]+dt_ext2:I_win_ext[-1]-dt_ext2]) 
		k=k+1 

	if in_['stack_norm'] : 
		D = s2d.norm_trans(np.array(cc['cc_mat']))
	else :
		D = np.array(cc['cc_mat'])
	(D,I1,I2)= s2d.stack_linear(D,in_['stack'],in_['stack_shift'])
	D        = s2d.filter(D,in_['p1'],in_['p2'],dt,npole=2)

	# spectral equalization done separatly on the causal/acausal side of the CC
	if in_['eqspec_ac'] : # seems less efficient
		(ndate, nwf) = D.shape 
		to = round(nwf/2)
		I_win_c = I_win[np.where(I_win > to)[0]]
		I_win_a = I_win[np.where(I_win < to)[0]]
		ref_fft_c = fft(ref[I_win_c])	
		ref_fft_a = fft(ref[I_win_a])	

		for idate in np.arange(0,ndate)  :
			Dc = D[idate,:]
			Dc_fft_c = fft(Dc[I_win_c])
			Dc_fft_a = fft(Dc[I_win_a])

			Dc_fft_c[np.where(abs(Dc_fft_c)==0)[0]]=0.001; 
			Dc_fft_a[np.where(abs(Dc_fft_a)==0)[0]]=0.001; 
			Dc_fft_c = Dc_fft_c / abs(Dc_fft_c) 
			Dc_fft_c = Dc_fft_c * abs(ref_fft_c)
			Dc_fft_a = Dc_fft_a / abs(Dc_fft_a) 
			Dc_fft_a = Dc_fft_a * abs(ref_fft_a)
			Dc[I_win_c] = np.real(ifft(Dc_fft_c))
			Dc[I_win_a] = np.real(ifft(Dc_fft_a))
			D[idate,I_win_c] = Dc[I_win_c]
			D[idate,I_win_a] = Dc[I_win_a]


	#equalization done on the full waveforms, i.e w/o windowing coda waves
	if in_['eqspec_full'] : 
		(ndate, nwf) = D.shape 
		to = round(nwf/2)
		#ref = ref / max(abs(ref))
		ref_fft = fft(ref)

		for idate in np.arange(0,ndate)  :
			Dc = D[idate,:]
			#Dc = Dc /max(abs(Dc))
			Dc_fft = fft(Dc)

			Dc_fft[np.where(abs(Dc_fft)==0)[0]]=0.001; 
			Dc_fft = Dc_fft / abs(Dc_fft) 
			Dc_fft = Dc_fft * abs(ref_fft)
			Dc = np.real(ifft(Dc_fft))
			D[idate,:] = Dc


	# classical spectral equalization performed on the coda window used to compute the dv/v
	if in_['eqspec'] : 
		(ndate, nwf) = D.shape 
		to = round(nwf/2)
		ref_fft = fft(ref[I_win])
		for idate in np.arange(0,ndate)  :
			Dc = D[idate,:]
			Dc_fft = fft(Dc[I_win])

			Dc_fft[np.where(abs(Dc_fft)==0)[0]]=0.001; 
			Dc_fft = Dc_fft / abs(Dc_fft) 
			Dc_fft = Dc_fft * abs(ref_fft)
			Dc[I_win] = np.real(ifft(Dc_fft))
			D[idate,I_win] = Dc[I_win]

			#if idate == 1 :
			#	import matplotlib.pyplot as plt
			#	#plt.plot(abs(ref_fft))
			#	#plt.plot(abs(Dc_fft),'r',linewidth=0.5)

			#	plt.plot(ref[I_win]/max(abs(ref[I_win])),color='black',linewidth=0.5)
			#	plt.plot(Dc[I_win]/max(abs(Dc[I_win])),color='red',linewidth=0.5)
			#	plt.plot(D[idate,I_win]/max(abs(D[idate,I_win])),color='blue',linewidth=0.5)
			#	plt.show()
			#	print(np.corrcoef(ref[I_win],Dc[I_win])[0,1])
			#	print(np.corrcoef(ref[I_win],D[idate,I_win])[0,1])
			#	ipdb.set_trace()




	
	#prepare output matrices and data : 
	[nstack,nwf]=np.shape(D) #D has changed this it has been stacked 
	out={} 
	out['dvv']=np.zeros((nstack))
	out['corrcoef']   =np.zeros((nstack))
	md_c={} 
	md_c['date1']= cc['date1'][I1]
	md_c['date2']= cc['date2'][I2]
	md_c['date'] = np.mean((md_c['date1'],md_c['date2']),axis=0) #to be checked 
	coeffmatrix=vcorrcoef2(refs[:,I_win],D[:,I_win]).transpose()

	# for convenience we extract directly the dvv for this path and the corres. corr coeff : 
	out['corrcoef'] =np.max(coeffmatrix,axis=1)
	#out['dvv']      =in_['stretch'][np.argmax(out['coeffmatrix'],axis=1)]
	out['dvv']      =in_['stretch'][np.argmax(coeffmatrix,axis=1)]

	if in_['print_dvv'] :
		print(out['dvv'])

	if in_['keep_coeffmatrix'] :
		out['coeffmatrix'] = coeffmatrix 

	if in_['export_matlab']  : 
		out_dir = 'export_matlab/dvv_'+str(in_['p1'])+'_'+str(in_['p2'])+'s'
		out_dir = out_dir+'_stack_'+str(in_['stack'])+'_'+str(in_['stack_shift'])
		#out_dir = 'dvv_matlab'
		if os.path.isdir(out_dir) == False :
			os.makedirs(out_dir)
		out_file = out_dir+'/'+cc['id'][0].decode()+'_'+cc['id'][1].decode()+'_'+str(int(round(cc['dist'])))+'km.mat'
		mat={}
		mat['dvv'] = out['dvv']
		mat['corrcoef']=out['corrcoef']
		D[:,0:200]=0
		D[:,-200:-1]=0
		D=np.transpose(s2d.norm(np.transpose(D)))
		mat['cc']=D#p.transpose(s2d.norm(np.transpose(D)))
		mat['ref']=ref 
		mat['t']=cc['t']
		mat['lon']=cc['lon']
		mat['lat']=cc['lat']
		mat['dist']=cc['dist']
		mat['I1'] = I_win_ext_start 
		mat['I2'] = int(win['neg']['coda'][1]) 
		mat['I3'] = int(win['pos']['coda'][0])
		mat['I4'] = I_win_ext_end 
		mat['date1']=md_c['date1']+366
		mat['date2']=md_c['date2']+366
		savemat(out_file,mat)


	return out, md_c, in_ 

def vcorrcoef(X,y):
    r_num = np.sum((X)*(y),axis=1)
    r_den = np.sqrt(np.sum((X)**2,axis=1)*np.sum((y)**2))
    r = r_num/r_den    
    return r 

def vcorrcoef2(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

    Returns
    -------
    np.array
      N X M array in which each element is a correlation coefficient.

    """

    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    #if n != y.shape[1]:
    #    raise ValueError('x and y must ' +
    #                     'have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    
    #cov = np.dot(x,y.T) - n * np.dot(mu_x[:, np.newaxis],mu_y[np.newaxis, :])
    cov = np.dot(x,y.T) # - n * np.dot(mu_x[:, np.newaxis],mu_y[np.newaxis, :])
    return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])


