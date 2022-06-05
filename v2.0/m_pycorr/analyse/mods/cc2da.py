from matplotlib.pyplot import * 
import matplotlib.patches as pt 
from matplotlib.collections import PatchCollection
#from bunch import Bunch as bunch 
import numpy as np
import pycorr_mods.dd as dd 
import m_pycorr.mods.s2d  as s2d
import pycorr_mods.lang as lang
try : 
	import ipdb
except :
	pass 

def get_window(trace,time_,dist,varargin={}) :
	""" INPUT : 
			- trace : single correlation trace 
			- time  : correlations time (lag) [s]
			- dist   : interstation distance  [km]
		OUTPUT : 
			- indice (not time!) of the beginning of the coda/sw in the positive/negative cc time

            NOTE : 
            Le cas des trajets longs n est pas parfaitement gere :
            On fait en sorte que la coda ne puisse pas aller au dela de nwf-5 
            par contre dans ce cas on verifie pas que le debut est avant nwf-5
            => peut induire des crash dans les codes suivants si nwf est trop petit 
            par rapport a la distance inter-station. 
	"""

	in_={}
	in_['sw_v1']  = 1.5   # Start of SW window [km/s]',...
	in_['sw_v2']  = 5.    # end of SW window [km/s];
	in_['hwidth'] = 2.    # half width of surface wave  [n*p2]
	in_['coda_dp']= 5.    # How many period after SW window, coda start [n*p2]    
	in_['p1']     = 5.    # corner period of the band pass filter [s]
	in_['p2']     = 10.   # idem 
	in_['plot']   =False  # plot window ?

	in_=lang.parse_options(in_,varargin)

	#get convenient variables here : 
	dt=time_[1]-time_[0]
	nwf=len(trace)
	hwidth_npts  = np.round(in_['hwidth']*in_['p2']/dt)
	coda_dp_npts = np.round(in_['coda_dp']*in_['p2']/dt)

	#filtering the trace between p1 and p2  : 
	trace=s2d.filter(trace,in_['p1'],in_['p2'],dt);

	#defining times where we look for surface waves
	t1=dist/in_['sw_v2'];
	t2=dist/in_['sw_v1'];
	if t2-t1 < dt    : # can happend for close stations
	  t2 = t1+5*dt     # subjective
	I0= np.where(time_ >=0)[0][0]
	#preparing output : 
	win    ={}
	win['pos']={}
	win['neg']={}
	win['success']=True 
	# si la distance trop grande => sw hors de la correlation :
	if t1 > time_[-1] :
		win['success']=False 
		return win 
	#analysing positive correlation time #on se laisse une marge de 5 pts :
	win['pos']['sw_guess']  = [np.nonzero(time_ >=t1)[0][0], np.nonzero(time_ <=t2)[0][-1] ]
	#for small distances : win['pos'][0] could equal to win['pos'][1] i.e dist/5 and dist1/.5 => point to the same indice
	win['pos']['sw_guess'][1] = max(win['pos']['sw_guess'][1],win['pos']['sw_guess'][0]+1)
	sw_I = np.argmax(abs(trace[win['pos']['sw_guess'][0]:win['pos']['sw_guess'][1]])) 
	win['pos']['sw_center']= sw_I+win['pos']['sw_guess'][0] 
	win['pos']['sw']=[]
	win['pos']['sw'].append(max((I0,win['pos']['sw_center'] - hwidth_npts)))
	win['pos']['sw'].append(min((nwf,win['pos']['sw_center']+hwidth_npts)))
	win['pos']['coda']=[]
	win['pos']['coda'].append(min((nwf-5,win['pos']['sw_center']+coda_dp_npts)))
	win['pos']['coda'].append(nwf-1)

	#analysing negative correlation time : 
	win['neg']['sw_guess']  = [np.nonzero(time_ >= -t2)[0][0],np.nonzero(time_ <= -t1)[0][-1]]
	#for small distances : win['neg'][0] could equal to win['neg'][1] i.e dist/5 and dist1/.5 => point to the same indice
	win['neg']['sw_guess'][1] = max(win['neg']['sw_guess'][1],win['neg']['sw_guess'][0]+1)
	sw_I=np.argmax(abs(trace[win['neg']['sw_guess'][0]:win['neg']['sw_guess'][1]]))
	win['neg']['sw_center']=sw_I + win['neg']['sw_guess'][0]
	win['neg']['sw']=[]
	win['neg']['sw'].append(max((1,win['neg']['sw_center'] - hwidth_npts)))
	win['neg']['sw'].append(min((I0,win['neg']['sw_center']+hwidth_npts)))
	win['neg']['coda']=[] 
	win['neg']['coda'].append(1) 
	win['neg']['coda'].append(max((5,win['neg']['sw_center']- coda_dp_npts)))

	if in_['plot'] == True : 
		fig,ax=subplots()
		plot(time_,s2d.norm(trace),color='k')
		for ica in ['pos','neg'] : 
			axvline(time_[win[ica]['sw_center']],color='red',lw=2)
			axvline(time_[win[ica]['sw'][0]]  ,  color='blue')
			axvline(time_[win[ica]['sw'][1]]  ,  color='blue')
			axvline(time_[win[ica]['coda'][0]],  color='yellow')
			axvline(time_[win[ica]['coda'][1]],  color='yellow')
			xlimit=np.round(dist*20)
			xlim((-xlimit,xlimit))
			ylim((-1.2,1.2))

			dx= time_[win[ica]['sw_guess'][1]] - time_[win[ica]['sw_guess'][0]]
			ax.add_patch(pt.Rectangle((time_[win[ica]['sw_guess'][0]],-1.2),dx,0.2,color='red',fill=True))

			dx= time_[win[ica]['sw'][1]] - time_[win[ica]['sw'][0]]
			ax.add_patch(pt.Rectangle((time_[win[ica]['sw'][0]],-1),dx,2,color='blue',fill=True,alpha=0.3))

			dx= time_[win[ica]['coda'][1]] - time_[win[ica]['coda'][0]]
			ax.add_patch(pt.Rectangle((time_[win[ica]['coda'][0]],-0.5),dx,1,color='yellow',fill=True,alpha=0.3))
		ipdb.set_trace()
		show()
	return win 
