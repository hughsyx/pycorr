import sys
import os 
import glob
import time
import warnings
import copy
import h5py as h5
#import silx.io as h5silx
import numpy as np 
import scipy.fftpack
import scipy.signal
import scipy.io as io
import skimage.restoration as deconv
from math import cos, sin, radians
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import offset_copy
import cartopy.crs as ccrs
from cartopy.io.img_tiles import Stamen, OSM
import cartopy.feature as cfeature
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics import locations2degrees
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size



def recursively_load_dic_from_h5(fin, path='/'):
    dic = {}
    for kname, item in fin[path].items():
        if isinstance(item, h5._hl.dataset.Dataset):
            dic[kname] = item[()]
            if dic[kname]==b'None':
                dic[kname] = None
        elif isinstance(item, h5._hl.group.Group):
            dic[kname] = recursively_load_dic_from_h5(fin, path + kname + '/')
    return dic



def mmap(sub,lon=[2.39,5.72],lat=[47.08,45.18],ev=[0.,0.],ch=None,proj='rob',filename='stations.png',info='station map',xylabel=False,res=None):
    if proj=='rob':
        ax = plt.subplot(sub,projection=ccrs.Robinson())
    elif proj=='azi':
        ax = plt.subplot(sub,projection=ccrs.LambertAzimuthalEqualArea(central_longitude=ev[0], central_latitude=ev[1]))
    else:
        ax = plt.subplot(sub,projection=ccrs.Robinson())            
    ax.set_global()
    ax.stock_img()
    ax.coastlines('110m')
    if ch == None:
        plt.scatter(lon, lat, c='red', s=20,transform=ccrs.PlateCarree(),cmap='plasma_r',alpha=0.5)
    else:
        lon = np.array(lon)
        lat = np.array(lat)
        all_ch_ind = ('111' == np.array(ch))
        Z_only     = ('100' == np.array(ch))
        others     = (1-all_ch_ind) * (1-Z_only)
        plt.scatter(lon[all_ch_ind], lat[all_ch_ind], c='red', s=20,transform=ccrs.PlateCarree(),cmap='plasma_r',alpha=0.5)
        plt.scatter(lon[Z_only], lat[Z_only], c='darkorange', s=20,transform=ccrs.PlateCarree(),cmap='plasma_r',alpha=0.5)
        plt.scatter(lon[others], lat[others], c='gold', s=20,transform=ccrs.PlateCarree(),cmap='plasma_r',alpha=0.5)    
    plt.title(info)
