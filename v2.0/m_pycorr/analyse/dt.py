################################################
# main_loop.py 
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
from scipy.fftpack     import fft
from scipy.fftpack     import ifft
#from scipy.signal      import correlate as xcorr 

try : 
	import ipdb
except :
	pass

from scipy.io import savemat 
	# tester linear ou spline ? 
	# modif de main loop pour retourner tt ce qu'il faut pour initializer un objet trace?

#ordre d'interpolation : 1 : insuffisant, 2 = OK, 3= super lent et pas d'interet ? 

def dt_ml(in_dir,in_dt,in_ml) : 
	in_={}
	in_['p1']=5              # used to filter + defining coda windows :
  	in_['p2']=10              
 	in_['v1']=1.5            # used to define coda windows :
  	in_['stack']=30          # how many date do we stack before measuring dvv ? 
  	in_['stack_shift']=1     # by how much do we shift the stack
 	#in_['limit_dt']=True     # look for dt less than p1/2 ?
 	in_['resampling']=100    # resampling factor 
 	in_['win_type'] ='sw'    # sw | coda | sw+coda | 
  	in_['coda_dp']=5         #  if coda : hm periods after sw, coda starts [n*p2]
  	in_['coda_duration']=120 #  if coda : coda end at in.tf 
  	in_['maxlag']       =10  # [s]
  	in_['sw_hwidth']   = 2.  # half width of surface wave  [n*p2]

	in_dt = lang.parse_options(in_,in_dt)
	tag='pydb_dt'+'__'+str(int(in_['p1'])).zfill(2)+'_'+str(int(in_['p2'])).zfill(2)+'s'
	tag=tag+'__stack_'+str(in_['stack'])+'_'+str(in_['stack_shift'])
	tag=tag+'__'+in_['win_type']
	if in_['win_type'] in ['sw'] : 
		tag=tag+'_'+str(in_['sw_hwidth'])
	if in_['win_type'] in ['coda','sw+coda'] :
		tag=tag+'__'+str(in_['coda_dp'])+'_'+str(in_['coda_duration'])+'s'
	in_['tag'] = tag 
	main_loop(in_dir,dt,in_dt,in_ml)


def dt(cc,in_dt) :
	in_={}
	in_['p1']=5              # used to filter + defining coda windows :
  	in_['p2']=10              
 	in_['v1']=1.5            # used to define coda windows :
  	in_['stack']=10          # how many date do we stack before measuring dvv ? 
  	in_['stack_shift']=1     # by how much do we shift the stack
 	#in_['limit_dt']=True     # look for dt less than p1/2 ?
 	in_['resampling']=100    # resampling factor 
 	in_['win_type'] ='sw'    # sw | coda | sw+coda | 
  	in_['coda_dp']=5         # if coda : hm periods after sw, coda starts [n*p2]
  	in_['coda_duration']=120 # if coda : coda end at in.tf 
  	in_['maxlag']       =10  # [s]
  	in_['sw_hwidth']   = 2.  # half width of surface wave  [n*p2]
	in_ = lang.parse_options(in_,in_dt)
#	ipdb.set_trace()

	#convenient variables : 
	ndate = len(cc['date1'])
	nwf   = len(cc['t'])
	dt    = cc['t'][1]-cc['t'][0]
	p2    = in_['p2']
	maxlag_npts = round(in_['maxlag']/dt)
	#determine the window used :
	ref=s2d.filter(cc['ref'],in_['p1'],in_['p2'],dt,npole=2);

	# apodize signal + window on selected wave 
	apo1  = get_cos2_window(nwf,in_['p2']*4,0,nwf,'inside')
	in_win={}
	in_win['p1']= in_['p1'] 
	in_win['p2']= in_['p2']
	in_win['plot']   = False 
	in_win['hwidth'] = in_['sw_hwidth'] 
	win   = cc2da.get_window(cc['ref'],cc['t'],cc['dist'],in_win)
	apo2  = get_window(win,in_['win_type'],nwf,p2)
	ref = ref *apo1*apo2 

	# linear stack of the daily CC : 
	(D,I1,I2)= s2d.stack_linear(cc['cc_mat'],in_['stack'],in_['stack_shift'])
	D        = s2d.filter(D,in_['p1'],in_['p2'],dt,npole=2)
	
	# useful variables to get the size of output matrices : 
	[nstack,nwf]=np.shape(D) #D has changed this it has been stacked 
	out={}
	out['coeff'] =  np.zeros((nstack))
	out['dt']    =  np.zeros((nstack))

	#main loop on stack : 
	#xcorr_time  = np.arange(-in_['maxlag'],in_['maxlag']+dt,dt)
	
	xcorr_time  = np.arange(-in_['maxlag'],in_['maxlag'],dt)
	xcorr_timei = np.arange(xcorr_time[0]+dt,xcorr_time[-1]-dt,dt/in_['resampling'])

	#dbg : 
	#mat_xcorr  = np.zeros((nstack,len(xcorr_time)))
	#mat_xcorri = np.zeros((nstack,len(xcorr_timei)))
	#mat_Dc     = np.zeros((nstack,nwf))
	for istack in np.arange(0,nstack) :
		D_current = D[istack,:] *apo1*apo2 
		if D_current.max() == 0 : 
			out['coeff'][istack] = 0
			out['dt'][istack]= 0 
			continue 
		xcorr = ctp_xcorr_norm(D_current,ref,maxlag_npts)
                #print str(len(xcorr_time))+' '+str(len(xcorr))
                #print str(maxlag_npts)
		print len(xcorr_time)
		print len(xcorr) 
		f_interp = interp1d(xcorr_time,xcorr,axis=0,kind=3)	
		xcorri=f_interp(xcorr_timei)
		I=np.argmax(xcorri)
		out['dt'][istack]=xcorr_timei[I]
		out['coeff'][istack]=xcorri[I]
		# dbg : 
		# mat_xcorr[istack,:]=xcorr 
		# mat_xcorri[istack,:]=xcorri 
		# mat_Dc[istack,:] = D_current 
	#dbg
	#mat={} 
	#mat['ref']=ref 
	#mat['D'] = D 
	#mat['dt']     =out['dt']
	#mat['coeff']  =out['coeff']
	#mat['ref_ori']=cc['ref']
	#mat['Dc']     =mat_Dc
	#mat['xcorr_time'] =xcorr_time 
	#mat['xcorr_timei']=xcorr_timei
	#mat['xcorr']  = mat_xcorr 
	#mat['xcorri'] = mat_xcorri
	#savemat('toto_'+str(round(cc['dist']))+'.mat',mat)

	
	md_c={} 
	md_c['date1']= cc['date1'][I1]
	md_c['date2']= cc['date2'][I2]
	md_c['date'] = np.mean((md_c['date1'],md_c['date2']),axis=0) #to be     
	return out, md_c, in_ 

	# main loop on days  

def get_window(win,win_type,nwf,p2) :
	#select correct portion of the signal + apodisation : 
	if win_type == 'coda'  :
		apo_pos=get_cos2_window(nwf,p2*2,win['pos']['coda'][0],win['pos']['coda'][1],'inside')
		apo_neg=get_cos2_window(nwf,p2*2,win['neg']['coda'][0],win['neg']['coda'][1],'inside')
	elif win_type == 'sw' :
		apo_pos = get_cos2_window(nwf,p2,win['pos']['sw'][0],win['pos']['sw'][1],'outside')
		apo_neg = get_cos2_window(nwf,p2,win['neg']['sw'][0],win['neg']['sw'][1],'outside')
	elif win_type =='sw+coda' : 
		apo_pos1=get_cos2_window(nwf,p2*2,win['pos']['coda'][0],win['pos']['coda'][1],'inside')
		apo_neg1=get_cos2_window(nwf,p2*2,win['neg']['coda'][0],win['neg']['coda'][1],'inside')
		apo_pos2 = get_cos2_window(nwf,p2,win['pos']['sw'][0],win['pos']['sw'][1],'outside')
		apo_neg2 = get_cos2_window(nwf,p2,win['neg']['sw'][0],win['neg']['sw'][1],'outside')
		apo_pos = apo_pos1+apo_pos2  
		apo_neg = apo_neg1+apo_neg2 
	apo = apo_pos+apo_neg 
	apo[np.where(apo>=1)]=1
	return apo 



def ctp_xcorr_norm(trace01, trace02,maxlag): # cross-coherence Time normalized before Corr
    lentrace   = len(trace01)
    maxlag= int(maxlag)
    goodnumber = int(maxlag+lentrace)

    tr2 = np.zeros(goodnumber)
    tr2[0:lentrace] = trace02
    tr2[0:lentrace] /= np.sqrt(np.sum(tr2[0:lentrace]**2))
    tr2 = fft(tr2,overwrite_x=True)
    tr2.imag *= -1
    tr1 = np.zeros(goodnumber)
    tr1[maxlag:int(maxlag+lentrace)]= trace01
    tr1[maxlag:int(maxlag+lentrace)] /= np.sqrt(np.sum(tr1[maxlag:int(maxlag+lentrace)]**2))
    tr2 *= fft(tr1,overwrite_x=True)
    return (ifft(tr2,overwrite_x=True)[0:2*maxlag+1].real)



# ------------------------------------------------------------------------
#
#
#               UTIlITY FUNCTIONS
#
#
#--------------------------------------------------------------------------
def get_cos2_window(nwf,napo,I1,I2,type='outside') :
	napo = int(napo)
	I1   = int(I1)
	I2   = int(I2)
	nwf  = int(nwf)
	cos2_up = np.power(np.cos(np.linspace(np.pi/2,0,napo)),2)
	cos2_dn = np.power(np.cos(np.linspace(0,np.pi/2,napo)),2)
	if type == 'inside' : 
		apo = np.zeros((nwf)) 
		gate = np.ones((I2-I1-napo*2))
		apo[I1:I2]=np.concatenate((cos2_up,gate,cos2_dn))
	elif type == 'outside' :
		apo = np.zeros((int(nwf+napo+napo))) 
		gate = np.ones((I2-I1))
		apo[I1:I2+napo*2]=np.concatenate((cos2_up,gate,cos2_dn))
		apo = apo[napo:-napo]
	return apo
