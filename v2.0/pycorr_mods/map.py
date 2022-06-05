import matplotlib.pyplot as plt
from cartopy.io.img_tiles import Stamen 
from obspy.geodetics import locations2degrees
import cartopy.crs as ccrs
import pickle
import numpy as np
import pycorr_mods.dd as dd 
import ipdb


def plot_station_map_from_dict_of_stations(db_sta,color='r',topo=False) :
	''' plot the stations of a pycorr data_5.0hz/daily/CH/2016/db.pkl file  
		
		typical usage is : 
		1) map.plot_station_map_from_dict_of_stations('data_5.0hz/daily/CH/2016/db.pkl',color='r',topo=False)
		OR 
		2) or using an already open db.pkl file :
			with open(db_path, 'rb') as f:
				db = pickle.load(f)
			map.plot_station_map_from_dict_of_stations(db['sta'],color='r',topo=False)
	'''

	if type(db_sta) == str : 
		with open(db_sta, 'rb') as f:
			db_sta = pickle.load(f)
			db_sta = db_sta['sta']
	edge = get_edge_from_dict_of_stations(db_sta,dx=3, dy=3) 
	background(map_edge=edge,topo=topo,coastline=True,res=None,draw_labels=True) 
	stations_from_dict_of_stations(db_sta,c=color,s=100,marker='^',alpha=1,edgecolors='k',label='',text={}) 


def get_edge_from_dict_of_stations(db_sta,dx=7, dy=7) :
	''' return map edge coordinate from a dict db_sta which contains one key per station : 
	db_sta['FR.ISO.00']['lat|lon']
	and so forth 
	'''
	min_lon = 1000 
	max_lon = -1000 
	min_lat = 100
	max_lat = -1000

	for ista in db_sta.keys() :
		min_lon=min(db_sta[ista]['lon'],min_lon)
		max_lon=max(db_sta[ista]['lon'],max_lon)
		min_lat=min(db_sta[ista]['lat'],min_lat)
		max_lat=max(db_sta[ista]['lat'],max_lat)
	#margin x and y :
	mx = (max_lon - min_lon)/dx 
	my = (max_lat - min_lat)/dy

	return [min_lon-mx,max_lon+mx,min_lat-my,max_lat+my]



def background(map_edge=[-10, 40, 40, 55],topo=False,coastline=True,res=None,draw_labels=True) :
	'''	in_file  : map_edge : coordinante of the 4 map edges [minlon,maxlon,minlat,maxlat] 
		topo     : True/False/1/0 : plot topography or not 
		coastline: True/False/1/0 : plot coastline or not 
		res	     : resolution of the topography : typically 1-8 
	'''

	# init new axis :
	plt.clf()
	ax = plt.axes(projection=Stamen('terrain-background').crs)
	ax.set_extent(map_edge,crs=ccrs.Geodetic())	

	if coastline : 	add_coastline(ax,map_edge)
	if topo : 
		if res==None :
			res = _get_topography_resolution(map_edge) 
			#res = get_resolution(map_edge)
			ax.add_image(Stamen('terrain-background'),res)
	gl = ax.gridlines(draw_labels=draw_labels) #xlabels_top=False,ylabels_right=False)
	gl.xlabels_top = False
	gl.ylabels_right = False
	plt.show()

def stations_from_dict_of_stations(sta,c='r',s=100,marker='^',alpha=1,edgecolors='k',label='',text={}) :
	''' to plot the station name, set a dictonnary with the following field for instance : 
		text = {} ;
		text['fs']=8 
		text['dy']=0.15 
		text['bbox'] = dict(boxstyle="round", ec=(1,1,1), fc=(1., 1, 1), alpha=0)
	'''
	ax = plt.gca()
	nsta = len(sta.keys())
	lon = np.zeros((nsta))
	lat = np.zeros((nsta))

	for k, ista in enumerate(sta) :
		lon[k] = sta[ista]['lon']
		lat[k] = sta[ista]['lat']

	if len(text) > 0 :
		for ista in sta.values() : 
			ax.text(ista['lon'],ista['lat']-text['dy'],ista['name'],fontsize=text['fontsize'],fontweight='bold',transform=ccrs.Geodetic(),ha='center',bbox=text['bbox']) 
	
	sc=ax.scatter(lon,lat,c=c,s=s,transform=ccrs.Geodetic(), marker=marker,alpha=alpha,edgecolors=edgecolors,label=label)



def add_coastline(ax,map_edge) :
	''' add a beautiful coastline to the map with a resolution that depends on map_edge '''
	minlon = map_edge[0]
	maxlon = map_edge[1]
	minlat = map_edge[2]
	maxlat = map_edge[3]
	r = int(locations2degrees(minlat,minlon,maxlat,maxlon)*111)
	if r <= 1000 :
		ax.coastlines('10m')
	if r > 1000 and r <= 5000:
		ax.coastlines('50m')    
	if r > 5000 and r <= 10000: 
		ax.coastlines('110m')
	if r > 10000: 
		ax.coastlines('110m')



def _get_topography_resolution(map_edge) :
	minlon = map_edge[0]
	maxlon = map_edge[1]
	minlat = map_edge[2]
	maxlat = map_edge[3]
	r = int(locations2degrees(minlat,minlon,maxlat,maxlon)*111)

	if r <= 1000 :
		res=8
	if r > 1000 and r <= 5000:
		res = 6
	if r > 5000 and r <= 10000: 
		res = 4
	if r > 10000: 
		res = 2
	return res 


