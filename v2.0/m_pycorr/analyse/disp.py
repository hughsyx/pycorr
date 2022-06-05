################################################
# disp.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


# A ameliorer : 
# re-tester systematiquement plusieurs max locaux => puis voir si en prenant un autre max on minimise la distance aux autres points 
from obspy import UTCDateTime
import m_pycorr.mods.dd as dd 
import m_pycorr.mods.lang as lang 
from m_pycorr.live.c1_trace import c1_trace
import m_pycorr.analyse.main_loop as main_loop 

import os
import glob
import obspy.core as obs 

from scipy.interpolate import interp1d 
from scipy.signal import argrelextrema
import numpy as np
from scipy import signal
from scipy import fftpack

try : 
    import ipdb
    from scipy.io import savemat 
except :
	pass


def disp_sisp_ml(in_dir,in_disp,in_ml) : 
	in_ = {}
	in_['periods']= np.arange(0.2,5.1,0.1)
	in_['nper']   = 100 
	in_['alpha']  = 0.08 
	in_['alpha2'] = 0.08 
	in_['alpha_dist']=[]
	in_['snr_p1'] = [0.2, 1, 3] 
	in_['snr_p2'] = [1  ,3, 5] 
	in_['snr_v1'] = 1 
	in_['snr_v2'] = 5
	in_['vmin']   = 1 
	in_['vmax']   = 5
	in_['use_local_max'] = True
	in_['tag']    = []
	in_['test_alpha'] = []
	#in_['test_alpha_dir']= 'test_disp_alpha'
	#in_['save_matrix']= False 
	in_disp       = lang.parse_options(in_,in_disp)

	in_ml['load_ref'] = True  
	in_ml['load_cc']  = False   
	# main loop options : for each trace read only the reference trace
	in_ml['load_ref'] = True  
	in_ml['load_cc']  = False   
	# defining output file name : in_['tag'] = dir name in_dir = full path
	tag = in_['tag']+'_'+	str(len(in_disp['periods']))+'_periods_a_'+str(in_['alpha'][0])
	if in_['use_local_max'] : 
		tag  = tag +'_ulm'
	in_disp['tag']   = 'pydb_sisp_'+tag 
	in_disp['in_dir']=  in_dir + '/'+in_disp['tag'] 
	# if test mode : define the output dir : 
	if len(in_disp['test_alpha']) > 0 :
		in_disp['test_alpha_dir'] = 'test_alpha/'+in_dir+'_'+in_disp['tag']
	main_loop.main_loop(in_dir,disp_sisp,in_disp,in_ml)

def disp_sisp(cc,in_disp) : 
	#sens est toujours NPS !!
	in_ = {}
	in_['periods']= np.arange(0.2,5.1,0.1)
	in_['nper']   = 100 
	in_['alpha']  = 0.08 
	in_['alpha2'] = 0.08
	in_['alpha_dist']=[]
	in_['snr_p1'] = [0.2, 1, 3] 
	in_['snr_p2'] = [1   ,3, 5] 
	in_['snr_v1'] = 1 
	in_['snr_v2'] = 5
	in_['vmin']   = 1
	in_['vmax']   = 5
	in_['use_local_max'] = True
	in_['tag']    = []
	#in_['save_matrix']=False 
	in_['test_alpha'] = []
	in_['test_alpha_dir']= '' 
	in_ = lang.parse_options(in_,in_disp)

	# create test_alpha_dir if necesssary 
	if len(in_['test_alpha']) > 0  :
		try : 
			os.makedirs(in_['test_alpha_dir'])
		except :
			pass 
	# first check that vg_fta is here : 
	if not os.path.isfile('vg_fta') : 
		dd.dispc('vg_fta cannot be found !! ','r','r')
		return
	# init common output : 
	md_c={}
	md_c['periods']= in_['periods']
	md_c['snr_p1']= in_['snr_p1']
	md_c['snr_p2']= in_['snr_p2']
	# init path specific output : SNR +GV 
	nper = len(in_['periods'])
	nper_snr = len(in_['snr_p1'])
	out={}
	out['gv']  = np.zeros((3,nper))     # velocity re-interpolated at discrete periods N/P/S
	out['snr'] = np.zeros((3,nper_snr)) # snr at each periods caus/acaus/sym
	
	# tmp dir specific to this process :
	pid  = str(os.getpid())
	cdir = os.getcwd() 
	out_dir  = in_disp['in_dir']+'/process_'+str(pid)
	out_sac = 'data.sac'

	dd.dd(cc)
	dd.dispc('len cc[ref] = '+str(len(cc['ref'])),'c','b')
	if cc['dist'] == 0 or cc['ref'].max == 0 :
		dd.dispc(' dist = 0 km or empty trace => exiting','r','b')
		return out, md_c, in_

	if not os.path.isdir(out_dir) : 
		os.makedirs(out_dir)
	os.chdir(out_dir)

	# convinient variable : 
	tau = cc['t'][1]-cc['t'][0]

	# init c1_trace instance to be improved :) : 
	cc['time']=cc['t']
	cc['tr'] = cc['ref']
	cc['tau']= tau;
	cc['title']=''
	c1_tr = c1_trace(cc)
	st_sac = c1_tr.sac_split_ca()
	
	# choose alpha value 
	if len(in_['alpha_dist']) == 0 :
		alpha = in_['alpha'][0]
		alpha2= in_['alpha2'][0]
	else :
		I_alpha=np.where(cc['dist'] > np.array(in_['alpha_dist']))[0][-1]
		alpha = in_['alpha'][I_alpha]
		alpha2= in_['alpha2'][I_alpha]
		dd.dispc('   dist = '+str(cc['dist'])+' : alpha = '+ str(alpha),'c','n')
		dd.dispc('                                alpha2 = '+ str(alpha2),'c','n')
	kca = -1
	for ica in ['N','P','S'] :
		kca = kca+1
		#st[0].data = np.array(c1[ica]['trace'])
		st_sac[kca].write(out_sac,format='SAC')
		# calling vg_fta :
		fmin = str(1./in_['periods'][-1])
		fmax = str(1./in_['periods'][1])
		cmd = cdir+'/vg_fta ' 
		cmd = cmd + out_sac +' fmin='+fmin+ ' fmax='+fmax+' vgMin='+str(in_['vmin'])	
		cmd = cmd + ' vgMax='+str(in_['vmax'])+' bmin='+str(alpha)+' bmax='+str(alpha2)
		cmd = cmd + ' disp=cont out=mat diag=PV'+' nfreq='+str(in_['nper'])+' ampMin='+'0'
		os.system(cmd)
	
		# read the P,V matrix, extracting the disp curve and interpolating it : 	
		dd.dispc('reading results','c','b')
		try :
			matrix = disp_sisp_read_matrix() 
			out['gv'][kca,:] = disp_sisp_interp_disp(matrix,in_['periods'],in_['vmin'],in_['vmax'],in_['use_local_max'])
			dd.dispc('done ','c','n')
		except : 
			dd.dispc('failed','r','r')

	# calcul des SNR : 	
	dd.dispc('computing SNR ','c','b')
	for iper in np.arange(0,nper_snr) : 
		c1_tr = c1_trace(cc) 
		#voir pk c'est necessaire dist trop faible => pb de fenetrage ?... 
		try : 
			in_snr={}
			in_snr['sw_v1']=in_['snr_v1']
			in_snr['sw_v2']=in_['snr_v2']
			snr=c1_tr.get_snr(in_['snr_p1'][iper],in_['snr_p2'][iper],in_snr)
			out['snr'][0,iper] = snr['neg']
			out['snr'][1,iper] = snr['pos']
			out['snr'][2,iper] = snr['sym']
			dd.dispc(str(iper)+' ok','g','r')
		except : 
			dd.dispc(str(iper)+' failed','r','r')
			pass 
	dd.dispc('done','c','n')

	# test_alpha : done on the symetric part of the CC 
	nalpha = len(in_['test_alpha'])
	if nalpha > 0 :
		kalpha = -1 
		for ialpha in in_['test_alpha'] :
			kalpha = kalpha + 1
			cmd = cdir+'/vg_fta ' 
			cmd = cmd + out_sac +' fmin='+fmin+ ' fmax='+fmax+' vgMin='+str(in_['vmin'])	
			cmd = cmd + ' vgMax='+str(in_['vmax'])+' bmin='+str(ialpha)+' bmax='+str(ialpha)
			cmd = cmd + ' disp=cont out=mat diag=PV disp=all'+' nfreq='+str(in_['nper'])+' ampMin='+'0'
			#PV : Period/Velocity
                        #FT : Frequency/Time
                        #PT : Period/Time
                        #FV : Frequency/Velocity
			os.system(cmd)

			matrix = disp_sisp_read_matrix() 
			gv     = disp_sisp_interp_disp(matrix,in_['periods'],in_['vmin'],in_['vmax'],in_['use_local_max'])

			if kalpha == 0 : 
				mnper = len(matrix['per']) 
				nperi = len(in_['periods'])
				mnvel = len(matrix['vel'])
				nalph = len(in_['test_alpha'])
				# init matlab dictionnary that will be saved : 
				mat={}
				mat['per']  = matrix['per']
				mat['vel']  = matrix['vel']
				mat['alpha']= in_['alpha'] 
				mat['alpha2']=in_['alpha2']
				mat['disp'] = np.zeros((mnper,mnvel,nalph))
				mat['gvi']  = np.zeros((nperi,nalph))
				mat['cc']   = cc
				mat['out']  = out 
				mat['in']   = in_
				#store result in the matlab dict : 
				dd.dispc('reading results','c','b')
				mat['disp'][:,:,kalpha] = matrix['disp']
				mat['gvi'][:,kalpha] = gv 
				dd.dispc('done ','c','n')
		#saving the result into a mat file : ]
		out_dir_alpha = cdir + '/'+in_['test_alpha_dir'] 
		out_mat= out_dir_alpha +'/'+str(np.round(cc['dist']*10)/10)+'_km_'+cc['id'][0]+'_'+cc['id'][1]+'.mat'
		if not os.path.isdir(out_dir_alpha) : 
			os.mkdir(out_dir_alpha) 
		savemat(out_mat,mat)
	os.system('rm *txt')
	os.system('rm *sac')
	os.chdir(cdir);
	return out, md_c, in_




def disp_sisp_read_matrix() : 
	matrix={} 
	matrix['disp'] = np.loadtxt('write_amp.txt')# freq, amp 
	matrix['vel'] = np.loadtxt('write_TV.txt') # velocity vector 
	matrix['per'] = np.loadtxt('write_FP.txt') # freq/per
	return matrix 

def  disp_sisp_interp_disp(matrix, periods,vmin,vmax,use_local_max) :
	try :
            if use_local_max  :
                #gv2 = matrix['vel'][matrix['disp'].argmax(axis=1)]
                gv = extract_gv_from_matrix(matrix,vmin,vmax)
            else :
                gv = matrix['vel'][matrix['disp'].argmax(axis=1)]

            I_per = np.where((periods >= matrix['per'][0]) & (periods <= matrix['per'][-1]))
            f_interp = interp1d(matrix['per'],gv,axis=0,kind='linear')
            gvi = np.zeros(len(periods))
            gvi[I_per] = f_interp(periods[I_per])
            return gvi
	except :
            gvi = np.zeros(len(periods))
            return gvi 

def extract_gv_from_matrix(matrix,vmin,vmax) : 
    # input : matrix['disp'] matrix['vel'] matrix['per']
    # analyse matrix['disp'], i.e the energy in the P-V domain by periods 
    # for each periods take the biggest local maximum between vmin and vmax excluding vmin and vmax
    dd.dispc(' extract_gv_from_matrix with local minimum stuff','g','r')
    #savemat('../../toto.mat',matrix)
    nper = matrix['disp'].shape[0] 
    vel = matrix['vel'] 
    I1 = np.where(vel > vmin)[0][0]
    I2 = np.where(vel < vmax)[0][-1]
    gv=np.zeros(nper)
    for iper in range(0,nper) :
        local_max_I  =argrelextrema(matrix['disp'][iper], np.greater)[0]
        local_max_val=matrix['disp'][iper][local_max_I]
        if len(local_max_I) > 0  : # si trace non nulle/platte :
            # sort the local maximum for greatest to the lowest 
            I_sort       = np.argsort(local_max_val)
            I_sort = I_sort[::-1] 
            local_max_I  = local_max_I[I_sort]
            local_max_val= local_max_val[I_sort]
            # look for the first local max which is not at I_sw1 nor I_sw2 
            xx=np.where((local_max_I > I1) & (local_max_I < I2))
            if len(xx) > 0 :
                if len(xx[0]) > 0  :
                    I_sw = xx[0][0] 
                    I_sw = local_max_I[xx[0][0]]
                else : #enveloppe creuse=> pas de max local hors bords  
                    I_sw = np.argmax(matrix['disp'][iper])
            else :
                I_sw = np.argmax(matrix['disp'][iper])
            gv[iper] = vel[I_sw] 
            #print gv[iper]
            #if gv[iper] < 0.6 :
            #    ipdb.set_trace()

    return gv 

def disp_sisp_read_disp(in_) : 
	disp=np.loadtxt('write_disp.txt')
	# reading the result : 
	try :
		per= disp[:,0] 
		gv = disp[:,1] 
	except :
		return False, False, False , False , False 
	#	pass 
	#	os.chdir(cdir);
	#	return out, md_c, in_

	I=np.argsort(per)
	per=per[I]
	gv =gv[I]
	# interpolating the periods : 
	I_per=np.where((in_['periods'] >= per[0]) & (in_['periods'] <= per[-1]))
	f_interp = interp1d(per,gv,axis=0,kind='linear')
	#out['gv'][kca,I_per] = f_interp(in_['periods'][I_per])
	gvi = f_interp(in_['periods'][I_per])
	return gvi, I_per, True, gv, per

#===========================================================================
#
#                  STANDARD YANG-STYLE DISP
#
#===========================================================================
def disp_ref_ml(in_dir,in_disp,in_ml) : 
	# dispersion curve measurements input : 
	in_           = {}
	in_['periods']= [5,10,20,40] 
	in_['save_matrix']=False 
	in_['alpha'] = False 
	in_['interp']= False # freq to which we interpolate [hz]
	in_disp       = lang.parse_options(in_,in_disp)
	# main loop options : for each trace read only the reference trace
	in_ml['load_ref'] = True  
	in_ml['load_cc']  = False   
	# defining output filename :
	tag='pydb_disp'+'__'+str(len(in_disp['periods']))+'_periods'
	if in_disp['save_matrix'] : 
		tag=tag+'_with_matrix'
	if in_['alpha'] : 
		tag = tag + '_alpha_'+str(in_['alpha'])
	in_disp['tag'] = tag 
	# lets go :
	main_loop.main_loop(in_dir,disp_ref,in_disp,in_ml)


def disp_ref(cc,in_disp) :
	in_={}
	in_['periods']=[5,10,20,40]
	in_['save_matrix']=False 
	in_['alpha']=False
	in_['interp']= False # [hz]
	in_ = lang.parse_options(in_,in_disp)

	#INTERPOLER ICI !
	if in_['interp'] != False :
		dd.dispc('  interpolating to '+str(in_['interp'])+' hz','c','n')
		f_interp = interp1d(cc['t'],cc['ref'],axis=0,kind='cubic')
		timei = np.arange(cc['t'][0],cc['t'][-1]+ 0.5/in_['interp'],1./in_['interp'])
		refi = f_interp(timei) 
		#cc['refo'] = cc['ref'].copy()
		#cc['to']= cc['t'].copy()
		cc['ref'] = refi
		cc['t']   = timei 
		#savemat('cc.mat',cc)

	fe = 1./(cc['t'][2]-cc['t'][1])
	
	out, env_matrix = measure_disp2(cc['ref'],cc['t'],cc['dist'],fe,1.0,5.0,in_['periods'],in_['alpha'])	

	if in_['save_matrix'] : 
		out['env_matrix'] = env_matrix 

	md_c={}
	md_c['periods']= in_['periods']
	return out, md_c, in_ 




def measure_disp2(cc,time,dist,fs,Vmin,Vmax,periods,alpha):
	''' measure dispersion curve + snr on a trace starting at the time 0
		trace : trace starting at the time 0 
		distance : (virtual) source-receiver distance in m
		f_sampling : sampling rate 
		Vmin, Vmax : velocity used to window the SW (km/s)
		periods  : list of periods at which the measurements is done

		TODO : voir SNR : 300 pts max ?
			 : gestion des faibles distances : lorsque I_vmin=I_vmax
	'''
	nper = len(periods)
	ntr  = len(cc)
	out={}
	out['snr']  = np.zeros((3,nper)) # snr at each periods acausal/causal/sym (?)
	out['tsw']  = np.zeros((3,nper)) # arrival time of the surface waves 
	out['peri'] = np.zeros((3,nper)) # instanteneous periods 
	out['gvi'] = np.zeros((3,nper)) # velocity measured at the instantenous fre 
	out['gv']  = np.zeros((3,nper)) # velocity re-interpolated at discrete periods 
	#out['per']  = periods            # out['per']

	# adjust the filter width to the distance 
	if  alpha == False :
		if (dist <200):   parameter = np.ones(150)*1.0    
		if (dist <500):   parameter = np.ones(150)*0.6   
		if (dist <1000):  parameter = np.ones(150)*0.4
		if (dist >=1000): parameter = np.ones(150)*0.3
	else :
		parameter = np.ones(alpha)*1.0 

    # compute the fft of the CC : 
	cc_fft = fftpack.fft(np.float32(cc))
	nper = len(periods) 
	freq = np.linspace(0.,fs,ntr)
	
	# defining temporal window to isolate SW : 
	I_sw1 = int(dist/Vmax*fs)
	I_sw2 = int(dist/Vmin*fs)
	I_sw2 = np.max([I_sw1+1,I_sw2])
	
	# number of points to be kept if storing the env : 
	env_matrix=np.zeros((3,nper,int(np.round(ntr/2))))
	#main loop on frequency : 
	for ii in range(0,nper) : 
		# apply a gaussian filter in the frequency domain  + symmetric correlation: 
		filter_gaussian = np.exp(-((freq-1./periods[ii])/(parameter[ii]*1./periods[ii]))**2)
		cc_fft_filtered = cc_fft*filter_gaussian
		cc_filtered = fftpack.ifft(cc_fft_filtered).real	
		# get symetric signal :
		cc_filtered_sym = np.fliplr([cc_filtered])+cc_filtered
		cc_filtered_sym = cc_filtered_sym[0]
		# get envelope and instantaneous frequency of the cc and sym cc : 
		cc_env, cc_freq        = get_env_instantenaous_freq(cc_filtered,fs)
		cc_sym_env, cc_sym_freq= get_env_instantenaous_freq(cc_filtered_sym,fs)
		
		#split causal/acausal enveloppe and correlation and add symmetric part :
		env_ca = split_ca(cc_env,time)
		env_ca['S']['trace']=split_ca(cc_sym_env,time)['P']['trace']
		#..
		cc_ca  = split_ca(cc_filtered,time) 
		cc_ca['S']['trace']= split_ca(cc_filtered_sym,time)['P']['trace']
		#.. 
		freqi_ca = split_ca(cc_freq,time)
		freqi_ca['S']['trace']= split_ca(cc_sym_freq,time)['P']['trace']

		#optional storage of the envelope : 
		env_matrix[0,ii,:]=env_ca['N']['trace']
		env_matrix[1,ii,:]=env_ca['P']['trace']
		env_matrix[2,ii,:]=env_ca['S']['trace']
		
		#beginning of the noise window :
		I_noise = int(np.min([ntr/2,I_sw2 + 2*periods[ii]]))
		#measure gv and snr on the pos/neg/enveloppe :
		kca =0 
		for ca in ['N','P','S'] :	
			env_ca[ca]['trace'][0:I_sw1]=0.
			env_ca[ca]['trace'][I_sw2:]=0.
			I_sw = np.argmax(env_ca[ca]['trace'])
			if (I_sw == I_sw1-1) or (I_sw== I_sw2-1) : #verifier indice ici
				# we look for the highest local maximum which is not exactly at the beg/end of the window
				local_max_I  =argrelextrema(env_ca[ca]['trace'], np.greater)[0]
				local_max_val=env_ca[ca]['trace'][local_max_I]
				# check that we have more than one local maximum
				if len(local_max_I) > 1  : 
					# sort the local maximum for greatest to the lowest 
					I_sort       = np.argsort(local_max_val)
					I_sort = I_sort[::-1] 
					local_max_I  = local_max_I[I_sort]
					local_max_val= local_max_val[I_sort]
					# look for the first local max which is not at I_sw1 nor I_sw2 
					try :
						I_sw = np.where((local_max_I > I_sw1) & (local_max_I < I_sw2-1))[0][0]
						I_sw = local_max_I[I_sw]
					except : #enveloppe creuse=> pas de max local hors bords  
						I_sw = np.argmax(env_ca[ca]['trace'])
			#else :
			#	if ca=='N' :
			#		dd.dispc(str(periods[ii])+'s '+ca+':  I_sw1 :'+str(I_sw1/5.)+' I_sw2 :'+str(I_sw2/5.)+' argmax :'+str(np.argmax(env_ca[ca]['trace'])/5.)+' I_sw :'+str(I_sw/5.),'c','n')
			t_sw = env_ca['time'][I_sw]
			if t_sw !=0 :
				out['gvi'][kca,ii]=dist/t_sw 
				out['tsw'][kca,ii]=t_sw 
			else  :
				out['gvi'][kca,ii]= 0. 
			out['snr'][kca,ii] = np.max(env_ca[ca]['trace'])/np.std(cc_ca[ca]['trace'][I_noise:])
			I1 = I_sw - 1*periods[ii]*fs
			I2 = I_sw + 1*periods[ii]*fs
			I1 = int(np.max([0,I1]))
			I2 = int(np.min([I2,ntr/2]))
			out['peri'][kca,ii]= 1./np.median(freqi_ca[ca]['trace'][I1:I2])
			kca=kca+1

	out['gv'] = out['gvi'].copy()
	for ica in [0,1,2] :
		I_per=np.where((periods >= out['peri'][ica,0]) & (periods <= out['peri'][ica,-1]))
		f_interp = interp1d(np.squeeze(out['peri'][ica,:]),np.squeeze(out['gvi'][ica,:]),axis=0,kind='linear')
		out['gv'][ica,I_per] = f_interp(periods[I_per])

	return out, env_matrix



#===========================================================================
#
#                      COMMON FUNCTIONS 
#
#===========================================================================

def get_env_instantenaous_freq(cc,fs):
	cc_analytic = signal.hilbert(cc)
	env = abs(cc_analytic)
	instantaneous_phase = np.unwrap(np.angle(cc_analytic)) 
	instantaneous_frequency = (np.diff(instantaneous_phase) /(2.0*np.pi) * fs)
	instantaneous_frequency = np.append(instantaneous_frequency,0);  # to have the same number of points
	return env,  instantaneous_frequency

def split_ca(tr,time) :
	m0=int(np.floor(len(time)/2))
	c1={}
	c1['time']  = time[m0:-1] 
	c1['P']={}
	c1['P']['trace'] = tr[m0:-1]
	c1['N'] ={} 
	c1['N']['trace'] = np.flipud(tr[1:m0+1])
	c1['S']={}
	c1['S']['trace'] = c1['P']['trace'] + c1['N']['trace']
	return c1 
