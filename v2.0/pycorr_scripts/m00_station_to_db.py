# This script generates the data_5.0hz/daily/net/year/db.mat and db.pkl files using a station coordinate file :
# datacenter net station_name loc lat lon elev depth 
# RESIF      FR  ISO          00  44.x 06.y 1200 0 

from obspy.core import UTCDateTime
import scipy.io as io
import copy
import pickle


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

def write_db_file(out_dir,in_,sta,ev={}) :
    ''' output a db.mat file containing all input parameters of the code as well as the station list
    The db.mat file is in data_5hz/daily/EU/2010 directory. date are converted to matlab format
    '''
    db_file      = out_dir+'/db.mat' 
    db_file_lock = db_file+'.lock'   
    #output db_file in a matlab-friendly way : 
    db                     = dict() 
    db['in_']              = copy.deepcopy(in_)
    db['in_']['day1']      = float(db['in_']['day1'].toordinal()) + 366.
    db['in_']['day2']      = float(db['in_']['day2'].toordinal()) + 366.
    if db['in_']['password'] is None:db['in_']['password']='None'
    if db['in_']['user_id'] is None:db['in_']['user_id']='None'
    if db['in_']['ev']['time_after'] is None:db['in_']['ev']['time_after']='None'
    db['sta']              = dict()
    db['sta']['lon']       = []
    db['sta']['lat']       = []
    db['sta']['sta']       = []
    db['sta']['net']       = []
    db['sta']['dc']        = []
    db['sta']['loc']       = []
    db['sta']['elev']      = []
    db['sta']['depth']     = []
    if ev is not {}:
        db['ev']           = copy.deepcopy(ev)
    for ista in sta :
        db['sta']['sta'].append(ista) 
        db['sta']['net'].append(sta[ista]['net'])
        db['sta']['dc'].append(sta[ista]['dc'])
        db['sta']['lon'].append(float(sta[ista]['lon']))
        db['sta']['lat'].append(float(sta[ista]['lat']))
        db['sta']['loc'].append(sta[ista]['loc'])
        db['sta']['elev'].append(float(sta[ista]['elev']))
        db['sta']['depth'].append(float(sta[ista]['depth']))
    io.savemat(db_file,db)
    #output db file in a python-friendly way :
    db        = dict()
    db['in_'] = in_
    db['sta'] = sta
    if ev is not {}:
        db['ev']           = copy.deepcopy(ev)
    filename  = db_file[0:-4]+'.pkl'
    ff        = open(filename,'wb')
    pickle.dump(db,ff)
    ff.close()



# -------------------------------------

def define_in_dict():
    in_=dict()
    in_['ch']                  = ['BH','HH']            # channel that we will attempt to download by descending priority 
    in_['cmp']                 = ['Z']                  # list of component we will download ([Z,N,E]) 
    in_['tag']                 = 'EU'                   # subdir were the data will be saved 
    in_['day1']                = UTCDateTime(2019,1,1)  # date range
    in_['day2']                = UTCDateTime(2021,12,31)  # 
    in_['format']              = 'h5'                   # only .h5 is supported 
    in_['station_list']        = 'station.txt'          # text file containing the list of station to be downloaded
    in_['path_out']            = './'                   # output path
    in_['sds_path']            = '/media/resif/validated_seismic_data/'
    in_['sds_metapath']        = '/media/resif/metadata/stationXML/'
    in_['user_id']             = None                   # user id for fdsn client
    in_['password']            = None                   # password for fdsn client
    in_['pp']                  = dict()                 # preprocessing  options : 
    in_['pp']['remove_resp']   = True                   # remove instrumental response ? not even loaded if SDS
    in_['pp']['cut_len']       = 0                      # old param
    in_['pp']['tap_len']       = 15*60                  #  - taper lenght in seconds, dwld request : 2*tap_len + 86400s  
    in_['pp']['freq']          = 25                     #  - frequency. to which data are decimated [hz] 
    try :
        in_['pp']['f_prefilt']     = (0.008, input_user['pp']['freq'] / 2 - 0.05 * input_user['pp']['freq'])
    except :
         in_['pp']['f_prefilt']     = (0.008, in_['pp']['freq'] / 2 - 0.05 * in_['pp']['freq'])
    in_['pp']['glitch']        = True                   # should we apply glitch correction ?
    in_['qc']                  = dict()                 # quality control options :
    in_['qc']['max_nsegments'] = 8       # - max. number of segement in each mseed file = max number of gaps 
    in_['qc']['min_duration']  = 72000   # - minimum length of signals per days 72000s=20h [s]
    in_['qc']['timeout']       = 300     # - maximum time to retrieve and process a single day [s]
    in_['ev']                  = {} # always '{}'
    in_['ev']['only_event']    = False # default False
    in_['ev']['event_list']    = 'events.txt' # input station file from 00_get_events.py or other
    in_['ev']['time_before']   = 0 # [s] from source time
    in_['ev']['time_after']    = None # [s] from source time, if None, autoscale from magnitude
    return in_



in_ = define_in_dict()
write_db_file('data_5.0hz/daily/FR/2019',in_,read_station_list('stationlist.txt'),ev={}) 

