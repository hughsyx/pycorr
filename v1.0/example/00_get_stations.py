import m_pycorr.m00_get_stations as get_sta
from obspy.core import UTCDateTime
in_                   = {}

###############
# Standard
###############  
# Here you'll find all possible options and their format 
# =====>    http://service.iris.edu/irisws/fedcatalog/1/
#in_['key']        = value 

in_['net']             = 'G' # 'net1,net2,net2...' accept wildcard
in_['sta']             = 'R*,T*' # 'sta1,sta2,sta3...' accept wildcard
in_['loc']             = '*' # 'loc1,loc2,loc3...' accept wildcard
in_['channel']         = 'LH*' # 'ch1,ch2,ch3...' accept wildcard
#in_['start']           = '2013-04-16' # yyyy-mm-dd
#in_['end']             = '2013-04-30' # yyyy-mm-dd
in_['startbefore']          = '2018-06-01' # yyyy-mm-dd
in_['endafter']             = '2018-08-31' # yyyy-mm-dd
#in_['lat']             = '45' # center point lat
#in_['lon']             = '0' # center point lon
#in_['maxradius']       = '10' # max radius
#in_['minradius']       = '0' # min radius
in_['includeavailability'] = True
# ...

stafile                = 'stations' # output name (.txt and .png)

############### RUN
get_sta.find_stations(input_user=in_,stafile=stafile)

