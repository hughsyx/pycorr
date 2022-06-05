################################################
# noise_count_station_per_day.py 
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 
import sys, h5py
import glob
import ipdb

def print_number_of_days(in_dir='data_1.0hz/daily/FR/2019', count_channel=False) :
	''' print in the terminal the number of stations per day :
		import noise 
		noise.show_number_of_days(net='CH',year=2016,in_dir='data_5.0hz') 
	'''
	# loop on days : 
	daily_files = glob.glob(in_dir+'/day*.h5')
	daily_files.sort()
	for iday in daily_files :
		nsta = 0 
		nchannel = 0 
		net_list=[]
		h5_day = h5py.File(iday,'r') 
		# now counting the number of channels : 
		for inet in h5_day : 
			if inet == '_metadata' : continue 
			nsta = nsta + len(h5_day[inet].keys())
			if count_channel :
				for ista in h5_day[inet] : 
					nchannel = nchannel +len(h5_day[inet][ista])
	
		if count_channel :
			dispc_help(iday+'  ',str(nsta).zfill(3)+' stations   ',str(nchannel).zfill(4)+' channels',30)#'c','n')
		else :
			dispc_help(iday+'  ',str(nsta).zfill(3)+' stations   ','',30)#'c','n')
		h5_day.close()
	

#------------------------------------------------
def dispc(text,color,attr) :
	# print text in color with attribute attr ! 
	if color =='r' :
		ncolor=31
	elif color=='g' :
		ncolor=32 
	elif color=='y' :
		ncolor=33
	elif color=='b' :
		ncolor=34
	elif color=='m' :
		ncolor=35
	elif color=='c' :
		ncolor=36
	elif color=='w' :
		ncolor=37
	elif color=='gray' :
		ncolor=38

	if attr=='b' :
		prefix='\x1b[1m';
	elif attr=='d' :
		prefix='\x1b[2m';
	elif attr=='u' :
		prefix ='\x1b[4m';
	elif attr=='blink' :
		prefix ='\x1b[5m'
	elif attr=='r' :
		prefix='\x1b[7m'
	else :
		prefix=''
	str_=prefix+'\x1b['+str(ncolor)+'m'+text+'\x1b[0m';
	print(str_)
    

def dispc_help(text1,text2,text3,length) :
    prefix=''#\x1b[1m';
    str_=prefix+'\x1b['+str(33)+'m'+text1.ljust(length,' ')+':'+'\x1b[0m';
    str_=str_+'\x1b['+str(37)+'m  '+text2.ljust(5,' ')+'\x1b[0m';
    str_=str_+'\x1b['+str(39)+'m'+text3[0:50]+'\x1b[0m';
    str_ = str_.replace('\n','')
    print(str_)



if __name__ =="__main__" :
	if len(sys.argv) == 2:
		print(sys.argv)
		print_number_of_days(in_dir=sys.argv[1]) 
	elif len(sys.argv)== 3 :
		print_number_of_days(in_dir=sys.argv[1],count_channel=True)
	else  :
		print('')
		print('USAGE : noise_count_station_per_day.py data_5.0hz/daily/CH/2016')
		print('')
	#noise.print_number_of_days(net='CH',year=2016,in_dir='data_5.0hz') 
