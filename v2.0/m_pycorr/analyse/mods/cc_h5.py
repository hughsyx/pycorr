import glob

def list(in_dir) :
	cc_list=glob.glob(in_dir+'/xcorr_*h5')
	cc_list.sort()
	return cc_list 

def get_path_number(ff) :
	return len(ff['md']['lon']) 
