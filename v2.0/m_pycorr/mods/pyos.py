################################################
# pyos.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


import os, glob 

def mkdir(dir_name) : 
	if not os.path.isdir(dir_name) : 
		try :
			os.makedirs(dir_name)
			return True 
		except :
			return False
	else :
		return False 
#---------------------
def remove(filename) : 
	try :
		os.remove(filename)
		return True 
	except :
		return False 

#----------------------------
def remove_dir(dir_name) : 
	file_list=glob.glob(dir_name+'/*')
	for ifile in file_list  :
		try :
			os.remove(ifile)
		except :
			return False
	try :
		os.removedirs(dir_name)
	except :
		return False 
	return True 
#-----------------------------------
def rename(src,dest) :
	try : 
		os.rename(src,dest)
		return True 
	except 	:
		return False 


#---------------------------------
def create_lock_file(filename) : 
	try :
		ff=open(filename,'wb')
		pickle.dump(filename,ff)
		ff.close()
		return True 
	except	:
		return False 
