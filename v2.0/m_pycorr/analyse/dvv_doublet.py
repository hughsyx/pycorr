import pycorr.mods.lang  as lang
import pycorr.mods.dd as dd
import pycorr.mods.s2d as s2d 
import pycorr.mods.pyos as pyos 
import pycorr.analyse.main_loop as main_loop 
import numpy as np 
from scipy.stats import linregress as linregress

from pycorr.live.c1_trace import c1_trace

# the following module are used for debugging/testing the code interactively: 
try :
	import ipdb 
	import matplotlib.pyplot as plt
	from scipy.io.matlab import savemat
	#plt.matshow(D)	
	#plt.show()
	#plt.plot([1,2,3,4])

except :	
	pass




# when comparing a current correlation with the reference, for each coda window : 
# - we split the coda into N windows of Y=in_['db_win'] [s] shifted by in_['db_shift] [s] 
# - we tab the coda window to  nfft points to get enough point in the freq domain. This is done w/o apodisation. 
#   apodising does not seem to improve the result. To be rechecked. 
# - the phase of the ref and current cc coda window are compared, assuming that the dt < to a half of period. 
#   we deal with the phase cycling problem by trying to add or sustract the phase by 2pi. We keep the phase difference #   which is the closest to zero. 
#
#  -We evaluate the correlation coefficient of each coda window between the ref and the current correlation. If it is 
#   greater than cc_min we use them for the linear-regression that establish dt= dvv*tau+b. This give a first 
#   measurement of dv/v
#  -We re-evaluate the dv/v using by computing for each coda window the corresponding dv/v_i and take the median value 
#   all  dvv/v_i . This make sense assuming there is no clock error in the data 
#  -we redo the same operation by removing coda window with cc > cc_min only and then again by further rejecting 
#   measurements > median dv/v +- one standard deviation.  
#
# To check the detail of the measurements i.e what happen to each coda window : 
# - set savemat_coda_window to True
# - then as it runs the coda generate some mat file containing all the information for a given couple of station, 
#   /date/coda window analysis. Results are stored in dvv_doublet/coda_window/
# - run pycorr.analysis.dvv_doublet_coda_win.m to plot the result in a fancy way ;-) 
#

def dvv_doublet(cc,varargin) : 
	in_={}
	in_['p1']= 1.             # measure dv/v between periods p1,p2 
	in_['p2']= 3.             #
	in_['v1']= 1.5            # minimum velocity of SW when looking for SW/coda 
	in_['coda_dp']= 5.        # coda starts after 5 periods followin the SW 
	in_['coda_length'] = 120. # assume coda has 120s of signel 
	in_['stack']= 10          # stack N correlations 
	in_['stack_shift']= 5.    # by which amount we shift the stacking window 
	in_['db_win'] = 30        # size of the coda window used for the dblt meas [s]
	in_['db_shift'] = 10      # by which amount we shift these windows [s]
	in_['cc_min']   = 0.8     # minimum correlation coefficient to keep a coda win  
	in_['nfft'] = 1024        # for each coda window is tabbed to get nfft point 
	in_['savemat_coda_window']  = False 
	in_['savemat_daily_meas']   = False 
	in_ = lang.parse_options(in_,varargin)

	def savemat_coda_window() : 
		mat={}
		mat['ref'] = ref 
		mat['D'] = D 
		mat['time'] = cc['time']
		mat['win1']=win1 
		mat['win2']=win2 
		mat['win_ref']=win_ref 
		mat['win_st'] =win_st
		mat['iwin']= iwin 
		mat['istack']=istack 
		mat['dphi'] = dphi 
		mat['freq'] = freq
		mat['dt']   = dt 
		mat['cc'] = meas['cc'][iwin]
		mat['cc_dt'] = meas['cc_dt'][iwin]
		mat['in'] = in_ 
		out_dir = 'dvv_doublet/coda_window' 
		pyos.mkdir(out_dir)
		filename = 'dist_'+str(int(ref.dist))+'km_'+ref.id[0]+'_'+ref.id[1]+'_win_'+str(iwin).zfill(3)+'.mat'
		savemat(out_dir+'/'+filename,mat)
		#print '.'

	def savemat_daily_meas() : 
		mat={}
		mat['ref'] = ref 
		mat['D'] = D 
		mat['time'] = cc['time']
		mat['win1'] = win1 
		mat['win2'] = win2 
		mat['istack']= istack 
		mat['meas']= meas 
		#mat['meas']= {}
		#mat['meas']['slope']= slope 
		#mat['meas']['intercept']=intercept 
		#mat['meas']['time'] = meas['time'][I_sel],
		#mat['meas']['dt']   = meas['dt_med'][I_sel]
		#mat['meas']['cc']   = meas['cc'][I_sel]
		#mat['meas']['std']  = meas['dt_std'][I_sel]
		mat['meas_all']={}
		mat['meas_all']['time'] = meas['time']
		mat['meas_all']['dt_med'] = meas['dt_med']
		mat['meas_all']['dt_mean'] = meas['dt_mean']
		mat['meas_all']['dt_std'] = meas['dt_std']
		mat['meas_all']['cc'] = meas['cc']

		mat['in'] = in_ 
		out_dir = 'dvv_doublet/daily_meas' 
		pyos.mkdir(out_dir)
		filename = 'dist_'+str(int(ref.dist))+'km_'+ref.id[0]+'_'+ref.id[1]+'_stack_'+str(istack).zfill(3)+'.mat'
		savemat(out_dir+'/'+filename,mat)

	# define convienient variable here : 
	ndate=len(cc['date1'])
	nwf  =len(cc['t'])
	tau = cc['t'][1]-cc['t'][0]
	
	# init the reference CC as a c1_trace instance ... to be improved :) : 
	cc['time']=cc['t']
	cc['tr'] = cc['ref']
	cc['tau']= tau;
	cc['title']=''
	ref = c1_trace(cc)
	ref.filter(in_['p1'],in_['p2'])

	# compute the coda window from the filtered reference :
	win_={} 
	win_['sw_v1'] = in_['v1'] 
	win_['sw_v2'] = 5  
	win_['coda_dp']     = in_['coda_dp']
	win_['coda_length'] = in_['coda_length']
	win = ref.get_window(in_['p2'],win_)

	# define the coda window where doublet method is applied separately on the pos/neg time : 
	win_pos1 = np.arange(win['pos']['coda'][0],win['pos']['coda'][1]-in_['db_win']+1,in_['db_shift'])
	win_pos2 = win_pos1  + in_['db_win']
	win_neg1 = np.arange(win['neg']['coda'][0],win['neg']['coda'][1]-in_['db_win']+1,in_['db_shift'])
	win_neg2 = win_neg1  + in_['db_win']
	win1 = np.concatenate((win_neg1,win_pos1))
	win2 = np.concatenate((win_neg2,win_pos2))

	# stack the daily correlations w/o filter 
	(D,I1,I2)= s2d.stack_linear(cc['cc_mat'],in_['stack'],in_['stack_shift'])
	D        = s2d.filter(D,in_['p1'],in_['p2'],tau,npole=2)

	#init the ouput structure 
	[nstack, nwf] = np.shape(D)
	out = {}
 	out['dvv_lg'] = np.zeros((nstack))   # slope of the linear-regression using selected window
 	out['dvv_lgb']= np.zeros((nstack))   # intercept [lock error] of the same linear regression 
 	out['dvv_lgr']= np.zeros((nstack))   # r_value of the linear regression 
 	out['dvv_lgp']= np.zeros((nstack))   # p_value of the linear regression 
	out['dvv_lgstd']= np.zeros((nstack)) # std_err of the linear regression 
 	out['dvv_med_all']     = np.zeros((nstack)) # median of the dt/t measured on all coda window 
 	out['dvv_med_all_std'] = np.zeros((nstack)) # the corresponding std 
 	out['dvv_med_sel']     = np.zeros((nstack)) # median dt/t on selected coda window
 	out['dvv_med_std']     = np.zeros((nstack)) # the corresponding std 
	out['dvv_med_sel_less_one_std'] = np.zeros((nstack)) # to the median dv/v +- one std 
#	out['dist']=np.zeros((nstack))+cc['dist'] # 
	#ipdb.set_trace()
	#out['id']  =cc['id']

	md_c={}
	md_c['date1']= cc['date1'][I1]
	md_c['date2']= cc['date2'][I2]
	md_c['date'] = np.mean((md_c['date1'],md_c['date2']),axis=0) #to be     
	md_c['p1']   = in_['p1'] 
	md_c['p2']   = in_['p2']
	md_c['stack']       = in_['stack']
	md_c['coda_length'] = in_['coda_length'] 
	md_c['stack']       = in_['stack']
	md_c['stack_shift'] = in_['stack_shift']
	md_c['db_win']      = in_['db_win'] 
	md_c['db_shift']    = in_['db_shift']

	# loop on the daily correlations. Meas = summary of the daily measurement performed on all coda-window
	nwin  = len(win1)
	for istack in np.arange(0,nstack): 
		if D[istack,:].max() == 0 : #no signal for this trace 
			continue 
		# we have data for this day : init a dict that will contain all coda win measurements :
		meas = {} 
		meas['dt_mean'] = np.zeros((nwin)) 
		meas['dt_med']  = np.zeros((nwin)) 
		meas['dt_std']  = np.zeros((nwin))
		meas['cc']      = np.zeros((nwin))
		meas['cc_dt']   = np.zeros((nwin))
		meas['time']    = np.zeros((nwin))
		# loop on coda windows & store the measurements for this window in meas
		for iwin in np.arange(0,nwin) :
			dt, dphi, win_ref, win_st, freq = \
				mesure_dt_this_coda_win(ref.tr, D[istack,:],tau,in_['nfft'],in_['p1'],in_['p2'],win1[iwin],win2[iwin]) 

 			meas['dt_mean'][iwin] = dt.mean()
 			meas['dt_med'][iwin]  = np.median(dt) 
 			meas['dt_std'][iwin]  = dt.std()
 			xcorr_coeff, xcorr_dt = vcorrcoef(ref.tr[win1[iwin]:win2[iwin]],D[istack,win1[iwin]:win2[iwin]])
 			meas['cc'][iwin]      = xcorr_coeff 
 			meas['cc_dt'][iwin]   = xcorr_dt*tau 
 			meas['time'][iwin]    = np.mean(ref.time[win1[iwin]:win2[iwin]])
 			
 			# optionnally save all the details of the measurements here : 
			if istack == 10 and in_['savemat_coda_window'] : 
				savemat_coda_window() 
			#else : 
 			#	dd.dispc('do not save ','g','b')
 		# Now intepret all coda-window measurement. Fist select coda windows we want : 
 		I_cc = np.argsort(meas['cc']) 
 		# on veut au moins les 5 meilleurs : 
 		I_bestcc= I_cc[np.arange(np.max((0,nwin-5)), nwin)]
 		# ou tous ceux avec un cc > 0.8 
 		I_cc = np.where(meas['cc'] > in_['cc_min']) #0.8) 
 		# apply both selection criteria and compute linear regression :
 		I_sel= np.sort(np.unique(np.concatenate((I_bestcc,I_cc[0]))))
 		slope, intercept, r_value, p_value, std_err = linregress(meas['time'][I_sel],meas['dt_med'][I_sel])
 		meas['I_sel'] = I_sel 
 		out['dvv_lg'][istack]    = slope 
 		out['dvv_lgb'][istack]   = intercept 
 		out['dvv_lgr'][istack]   = r_value 
 		out['dvv_lgp'][istack]   = p_value 
 		out['dvv_lgstd'][istack] = std_err 
 		# assuming there is no clock error : take the median dt/tau on all meas + selected
 		dvv_sel = meas['dt_std'][I_sel]/meas['time'][I_sel] 
 		dvv_std = np.std(dvv_sel)
 		dvv_med = np.median(dvv_sel)
 		out['dvv_med_all'][istack]     = np.median(meas['dt_std']/meas['time']) # median dt/t on all coda window 
 		out['dvv_med_all_std'][istack] = np.std(meas['dt_std']/meas['time'])    # the corresponding std 
 		out['dvv_med_sel'][istack] = dvv_med                                    # median dt/t on selected coda window
 		out['dvv_med_std'][istack] = dvv_std                                    # the corresponding std 
 		I=np.where( abs(dvv_sel - dvv_med) < dvv_std)[0]               # median dt/t on selected coda window
 		if len(I)> 0 :                                                 # reselected to keep only dt/t close to 
 			out['dvv_med_sel_less_one_std'][istack] = np.median(dvv_sel[I])   # to the median dv/v +- one std 
 		else :
 			out['dvv_med_sel_less_one_std'][istack] = dvv_med
 		# optionnally save the daily measurements into a matlab file : 
		if in_['savemat_daily_meas'] :
			savemat_daily_meas() 

	return out, md_c, in_

def mesure_dt_this_coda_win(ref,cc,tau,nfft,p1,p2,I1,I2) : 
	''' mesure delta between the reference trace ref, and the current correlation cc between the periods p1 and p2. 
		I1 and I2 are the indice of the beg/end of the coda window. This function return variable useful for doublet analysis, and to check the measurements (savemat) '''

	# first define the indice of the frequency where we measure the dv/v 
	win_npts  = I2-I1+1 
	freq = np.fft.rfftfreq(nfft)/tau
	I_f1 = np.where(freq >= 1./p2)[0][0]
	I_f2 = np.where(freq <= 1./p1)[0][-1]

	# win_ref = the reference trace  tabbed to nfft point (to have enough frequency points)
	win_ref=np.zeros((nfft))
	win_ref[0:win_npts-1] = ref[I1:I2] 
	ref_fft=np.fft.rfft(win_ref)
	
	# win_cc = the current correlation trace  tabbed to nfft point (to have enough frequency points)
	win_cc=np.zeros((nfft)) 
	win_cc[0:win_npts-1] = cc[I1:I2] # D[istack,db_pos1[iwin]:db_pos2[iwin]]
	cc_fft=np.fft.rfft(win_cc)

	# compare the phase of the ref +- 2pi and the current cc. 
	nphi = I_f2- I_f1
	dphi1 = np.zeros((3,nphi))
	dphi1[0,:] = np.angle(ref_fft[I_f1:I_f2]) - np.angle(cc_fft[I_f1:I_f2]) 
	dphi1[1,:] = dphi1[0,:] + 2*np.pi  
	dphi1[2,:] = dphi1[0,:] - 2*np.pi 
	dphi = np.zeros(nphi) 

	#keep the phase difference closes to zero and compute the corresponding dt
	for iphi in np.arange(0,nphi) :
		dphi[iphi] = dphi1[abs(dphi1[:,iphi]).argmin(),iphi]
		dt = dphi/(2*np.pi)#/freq
		dt = dt / freq[I_f1:I_f2]

	return dphi, dt , win_ref, win_cc, freq[I_f1:I_f2] 



	# --loop on the coda window 
	# --taper + phase difference + cc ?
	# -- get dt as a function of the time 
	# -- linear regression 
	# -- store ?
	#plt.matshow(D)

def vcorrcoef(X,y):
	xcorr = np.correlate(X,y,'same')
	r_dt  = xcorr.argmax()  - (len(xcorr)/2)
	r_num = xcorr.max() 
	r_den = np.sqrt(np.sum((X)**2,axis=0)*np.sum((y)**2))
	r = r_num/r_den    
	#r_num = np.sum((X)*(y),axis=0)a

	#r_den = np.sqrt(np.sum((X)**2,axis=0)*np.sum((y)**2))
	#r = r_num/r_den    

	return r, r_dt 



def dvv_ml(in_dir,in_dvv,in_ml):
	in_={}
	in_['p1']= 1.     
	in_['p2']= 3. 
	in_['v1']= 1.5 
	in_['coda_dp']= 5. 
	in_['coda_length'] = 120. 
	in_['stack']= 10 
	in_['stack_shift']= 5.
	in_['db_win'] = 30   # size of the coda window used for the dblt meas [s]
	in_['db_shift'] = 10 # by which amount we shift these windows [s]
	in_['cc_min']   = 0.8     # minimum correlation coefficient to keep a coda win  
	in_['nfft'] = 1024 #
	in_['savemat_coda_window']  = False 
	in_['savemat_daily_meas']   = False 

	in_dvv = lang.parse_options(in_,in_dvv)
	dd.dispc('computing dvv using the doublet method on '+in_dir,'y','b')
	dd.dd(in_dvv)

	tag = 'pydb_dvv_doublet_'+'__'+str(int(in_['p1'])).zfill(2)+'_'+str(int(in_['p2'])).zfill(2)+'s'
	tag=tag+'__stack_'+str(in_['stack'])+'__'+'stack_shift_'+str(in_['stack_shift'])
	tag=tag+'__coda_length_'+str(int(in_['coda_length']))
	in_dvv['tag']=tag;

	dd.dispc('file name will prefixed with : '+tag,'g','b')

	main_loop(in_dir,dvv_doublet,in_dvv,in_ml) 


	