# This function 
#  - reads the full C3 files (i.e containing the contribution of each Vs to the correlation) 
#  - for each pair of station, select virtual sources in the good range of azimuth and distance 
#  - average their contribution (i.e the C2  corresponding to each selected source). 
#  - output the result using the same file format as the C1 
#
# Could be reformalized and integrated inside the pycorr package directly.  

#general
import os , glob , h5py, pickle, random 
# pycorr 
#import pycorr_x.modules.h5 as h5 
import m_pycorr.mods.h5 as h5
import m_pycorr_c3.live.c3 as c3 
import m_pycorr.mods.dd as dd 
import m_pycorr.mods.lang as lang 
import numpy as np

try : 
	import ipdb
except :
	pass 

#Work with ZZ C3 only (when converting md_c from c3 to c1 + when calling h.c2_to_c1(I[ipath],0,in_vs,in_)

def c3_to_c1(in_dir,varargin_vs_sel={},varargin_c3_to_c1={}) : 
	# options that contol the way virtual sources are selected :
	in_vs= {} 
	in_vs['dist'] = [10, 1000] 
	in_vs['daz']  = 5 
	in_vs['single_source']=False
	in_vs['rm_btw_sta'] = True  
	in_vs=lang.parse_options(in_vs,varargin_vs_sel)
	# options that controls how they are summed etc :
	in_ = {}
	in_['savemat'] = False # only for testing purpose 
	in_['norm']    = True     # normalize the contribution of each vs b4 summing them ?
	in_['tag']     = '' 
	in_ = lang.parse_options(in_,varargin_c3_to_c1)

	# ex dict will contain all information about the flow of the code : 
	ex = {}      
	ex['in_']     = in_ 
	ex['in_vs']   = in_vs
	ex['out_dir'] = ex_get_output_dirname(ex)

	# now use the class c3_live to open the C3 and loop on the C3/xcorr*h5 files :

	h = c3.c3_live(in_dir)
	nfile = len(h.file) 
	c3_I_file = [*range(nfile)]
	print(nfile)
	print(c3_I_file)
	random.shuffle(c3_I_file)
	#c3_file_list = h.file 
	#random.shuffle(c3_file_list)
	for I_c3 in c3_I_file :
		ic3 = h.file[I_c3]
		ex_ic3 = ex_get_output_filename(ex,ic3) 
		# should we work on this xcorr file :
		if len(glob.glob(ex_ic3['h5_file_search'])) > 0 :
			dd.dispc('  '+ex_ic3['h5_file_search']+' exist => continuing','r','b')
			continue 
		if len(glob.glob(ex_ic3['lock_file_search'])) > 0 :
			dd.dispc('  '+ex_ic3['lock_file_search']+' exist => continuing','r','b')
			continue 
		if os.path.isfile(ex_ic3['out_file']):
			dd.dispc('  '+ex_ic3['out_file']+' exist => continuing','r','b')
			continue 

		# lets work in this xcorr file : 	
		ex_create_lock_file(ex_ic3['lock_file'])  
		dd.dispc('  working on '+ic3,'g','b')

		# add metadata and tweak them so that it looks a normal c1 file :
		ff_in  = h5py.File(ic3,'r')
		ff_out = h5py.File(ex_ic3['h5_file'],'w') 
		ff_out.copy(ff_in['md'],'/md')
		ff_out.copy(ff_in['md_c'],'/md_c')
		del ff_out['md_c']['c1_dir']
		del ff_out['md_c']['t_c1a']
		del ff_out['md_c']['t_c1b']
		ff_out['md_c']['date1_ordinal'] = ff_out['md_c']['date1']
		ff_out['md_c']['date2_ordinal'] = ff_out['md_c']['date2']
		ff_out['md_c']['cmp']='ZZ'
		
		h5.copy_dict_to_group(ff_out,{'in_':in_})
		ff_out.create_dataset('/in_/cc_cmp', data = 'ZZ')
		
		# loop on path : 
		path_list = range(int(h.file_I[I_c3,0]),int(h.file_I[I_c3,1])+1)
		npath_list = len(path_list)
		ff_out.create_dataset('/ref_nstack', data = np.ones((1,npath_list)))
		for ipath in path_list :
			mat = h.c2_to_c1(ipath,0,in_vs,in_)
			dset_root = '/'+str(h.id[ipath][0],'utf8')+'/'+str(h.id[ipath][1],'utf8')+'/ZZ'
			ff_out.create_dataset('/ref'+dset_root, data= mat['cc_sum'])
			
		ff_out.close()
		ff_in.close()
		# renaming everything 
		try :
			os.rename(ex_ic3['h5_file'],ex_ic3['out_file'])
		except :
			path 
		try :
			os.remove(ex_ic3['lock_file']) 
		except :
			pass 


def ex_create_lock_file(filename) : 
	a=4; 
	ff=open(filename,'wb')
	pickle.dump(a,ff)
	ff.close()

def ex_get_output_filename(ex,ic3) :
	xcorr_name = ic3.split('/')[1]
	root_name = xcorr_name[0:-3]
	pid = str(os.getpid())
	ex_ic3 = {} 
	ex_ic3['out_file']  = ex['out_dir']+'/'+xcorr_name
	ex_ic3['lock_file'] = ex['out_dir']+'/'+root_name+'_'+pid+'.lock'
	ex_ic3['h5_file']   = ex['out_dir']+'/'+root_name+'_'+pid+'.h5'
	ex_ic3['lock_file_search'] = ex['out_dir']+'/'+root_name+'_*.lock'
	ex_ic3['h5_file_search']   = ex['out_dir']+'/'+root_name+'_*h5'
	return ex_ic3

def ex_get_output_dirname(ex) : 
	if ex['in_vs']['single_source']==False :
		out_dir = 'C3_to_C1__dist_'+str(ex['in_vs']['dist'][0])+'_'+str(ex['in_vs']['dist'][1])+'km'
	else :
		out_dir = 'C3_to_C1__single_source__dist_'+str(ex['in_vs']['dist'][0])+'_'+str(ex['in_vs']['dist'][1])+'km'
	out_dir = out_dir +'__daz_'+str(ex['in_vs']['daz'])
	if ex['in_']['norm'] : 
		out_dir = out_dir +'__norm'
	if ex['in_vs']['rm_btw_sta'] :
		out_dir = out_dir +'__rm_btw_sta'
	out_dir = out_dir + ex['in_']['tag']
	mkdir(out_dir)
	return out_dir 

def mkdir(dir_name) :
	try :
		os.makedirs(dir_name)
	except :
		pass 


# restricted config : 
#c3_to_c1('C2')


in_vs= {} 
in_vs['dist'] = [5, 60] 
in_vs['daz']  = 5
in_vs['rm_btw_sta'] = True  
in_vs['single_source']=False

in_={}
in_['max_pos']=False
c3_to_c1('C3__m50_compute_c3_new_code_test_single_vs_01-05s__sw',in_vs,in_)
