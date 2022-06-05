################################################
# station_plot.py
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# Sept 2020
################################################


#!/usr/bin/env python3 

import pycorr_mods.map as map_
import pycorr_mods.dd as dd
import ipdb 


def plot_station_list(filename,dx=7,dy=7,topo=False) : 
	sta = read_station_list(filename)
	map_edge = map_.get_edge_from_dict_of_stations(sta,dx=dx,dy=dy) 
	ax=map_.background(map_edge=map_edge,topo=topo,coastline=True,res=None)
	map_.stations_from_dict_of_stations(sta)
	


def read_station_list(filename) :
    sta={}
    ff=open(filename,'r')
    for iline in ff : 
        (dc,net,name,loc,lat,lon,elev,depth) = iline.split()
        if loc=='--' : loc = ''
        if loc == u'' : kname=net+'_'+name+'_00'
        else: kname=net+'_'+name+'_'+loc 
        sta[kname]={}
        sta[kname]['net']  = net
        sta[kname]['name'] = name
        sta[kname]['loc']  = loc
        sta[kname]['kname']= kname
        sta[kname]['lat']  = float(lat)
        sta[kname]['lon']  = float(lon)
        sta[kname]['elev'] = float(elev)
        sta[kname]['depth']= float(depth)
        sta[kname]['dc']   = dc
    ff.close()
    return sta 



if __name__=="__main__" :
	print('x')
	if len(sys.argv)>2 :
		filename = sys.argv[1] 
		plot_station_list(filename,dx=7,dy=7,topo=False) 
	else :
		dd.dispc('USAGE : station_list_plot stations/2018_FR.txt topo=1','w','n')



