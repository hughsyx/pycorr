################################################
# m00_get_stations.py 
# Pierre Boue (UGA)
# Laurent Stehly (UGA)
# email: pycorr@univ-grenoble-alpes.fr
# May 2017
################################################
"""

"""

import pickle
import urllib
import os
import numpy as np 
import m_pycorr.mods.lang as lang
import m_pycorr.mods.dd as dd
from obspy.core import UTCDateTime
from obspy.geodetics import locations2degrees

global plotmap 
plotmap = False
try:
  import matplotlib.pyplot as plt
  from matplotlib.transforms import offset_copy
  import cartopy.crs as ccrs
  from cartopy.io.img_tiles import Stamen
  import cartopy.feature as cfeature
  plotmap = True
except:
  pass

#import m_pycorr.mods.dd as dd
import ipdb


##########################################################################################
##########################################################################################
##          ##             #####    #####             ##    ##          ##       ####   ##
#    #########             ####      ####             ##    ##          ##       ####   ##
#    ##############   ########        ########   #######    ##   ####   ##   ##   ###   ##
##       ##########   #######    ##    #######   #######    ##   ####   ##   ##   ###   ##
######      #######   ######    ####    ######   #######    ##   ####   ##   ###   ##   ##
##########   ######   #####              #####   #######    ##   ####   ##   ###   ##   ##
##########   ######   ####    ########    ####   #######    ##          ##   ####       ##
##          #######   ###    ##########    ###   #######    ##          ##   ####       ##
##########################################################################################
##########################################################################################


def find_stations(input_user_fdsn={},input_user={}):
    ''' this function check the data availability on FDSN database 
        It creates the *stafile*.txt file that is needed for following processing
    '''

    in_fdsn                   = {}
    in_fdsn['net']            = 'FR,GR'
    in_fdsn['sta']            = '*'
    in_fdsn['loc']            = '*'
    in_fdsn['channel']        = 'BH*'
    in_fdsn['format']         = 'text'
    in_fdsn['includeoverlaps'] = 'true' #if the same channel is available at two dc, output it. Should be True 
    in_fdsn['nodata']          = '404'
    in_fdsn = lang.merge_options(in_fdsn,input_user_fdsn)

    in_={}
    in_['lat']=[-90, 90] 
    in_['lon']=[-360, 360] 
    #in_['bb_only']=False # 
    in_['sp_but_no_bb'] = False # list stations is short periods channels [E,S,M]HZ are available but no broadband
    in_['acc_but_no_vel'] = True   # list stations if accelerometric channels [*N*] data are available but no vel
    in_ = lang.merge_options(in_,input_user)

    dd.dispc('Entering get station','y','b')
    print_documentation(in_fdsn,in_)

    url      = 'http://service.iris.edu/irisws/fedcatalog/1/query?'
    data     = urllib.parse.urlencode(in_fdsn)
    f    = urllib.request.urlopen(url + data)
    sta_info    = f.read().decode('utf-8').split('\n')



    sta_dict = sta_info_to_dict(sta_info,in_) 
    sta_dict = sta_dict_remove_empty_net_dc_sta(sta_dict)
    sta_dict_count(sta_dict,'initially we found :')

    #if requested by user keep only stations is there is a short period channel and no broadband 
    if in_['sp_but_no_bb'] == True :    
        sta_dict_filter_channel(sta_dict,in_) 
        sta_dict = sta_dict_remove_empty_net_dc_sta(sta_dict)
        sta_dict_count(sta_dict,'After keeping only SP if BB are not here :')

    #if requested by user keep only stations if there is an accelerometers and no velocimeter 
    if in_['acc_but_no_vel'] == True :    
        sta_dict_filter_channel(sta_dict,in_) 
        sta_dict = sta_dict_remove_empty_net_dc_sta(sta_dict)
        sta_dict_count(sta_dict,'After keeping only accelerometers if velocimeters are not here :')

    # output it into several files with different kind of organization :
    sta_dict_to_pycorr_per_year(sta_dict,in_fdsn)
    sta_dict_to_pycorr_per_year_and_network(sta_dict,in_fdsn) 
    sta_dict_to_pycorr_per_year_and_network_full_info(sta_dict,in_fdsn) 
    sta_dict_to_pycorr_bulk(sta_dict,in_fdsn) 

def sta_dict_to_pycorr_bulk(sta_dict,in_fdsn) :
    ''' here we just put all stations in one file with only one channel per station'''

    out_dir = 'stations_bulk'
    if not os.path.isdir(out_dir) :
        os.mkdir(out_dir)

    ff = open(out_dir+'/station_all.txt', "w") # create the output file 
    for inet in sta_dict :                                 # loop on the network/dc/station/ch 
        for idc in sta_dict[inet] : 
            for ista in sta_dict[inet][idc] :
                for ich in sta_dict[inet][idc][ista] :
                    csta = sta_dict[inet][idc][ista][ich]
                    cline = csta['dc']+'   '+csta['net']+'   '+csta['sta'] +'   '+csta['loc']+'   ' + csta['lat']
                    cline = cline + '   '+csta['lon']+'   '+csta['elev']+'   '+csta['depth']+'\n'
                    ff.write(cline)
                    break # we output only the station one time (and not for each channel)
    ff.close()
   

def sta_dict_to_pycorr_per_year_and_network_full_info(sta_dict,in_fdsn) :
    ''' here we sort the stations we found per year and network, i.e we create one file per year/network
        if a station has two datacenter, we output it twice 
        We put all stations informations !'''

    out_dir = 'stations_per_year_network_full_info'
    if not os.path.isdir(out_dir) :
        os.mkdir(out_dir)

    year1 = in_fdsn['start'].split('-')[0]
    year2 = in_fdsn['end'].split('-')[0]

    for iyear in np.arange(int(year1),int(year2)+1) : # loop on year 
        for inet in sta_dict :                                 # loop on the network/dc/station/ch 
            filename=out_dir+'/'+str(iyear)+'_'+inet+'.txt'
            ff = open(filename, "w") # create the output file 
            kline=0
            for idc in sta_dict[inet] : 
                for ista in sta_dict[inet][idc] :
                    for ich in sta_dict[inet][idc][ista] :
                        csta = sta_dict[inet][idc][ista][ich]
                        if int(csta['start'].split('-')[0]) > iyear : # si la stationa demarre apres l'anne courrante
                            continue 
                        if int(csta['end'].split('-')[0]) < iyear : # si la station s'est arrete avant l'anne courrante
                            continue 
                        cline = csta['dc']+'   '+csta['net']+'   '+csta['sta'] +'   '+csta['loc']+'   ' + csta['lat']
                        cline = cline + '   '+csta['lon']+'   '+csta['elev']+'   '+csta['depth']+'  '+csta['ch']
                        cline = cline +'   '+csta['start']+ '  '+csta['end']+'   '+csta['dip']+'   '+csta['azimuth']
                        cline = cline +'   '+csta['sensor'].replace(' ','_')+'\n'
                        ff.write(cline)
                        kline = kline+1
                        #reak # we output only the station one time (and not for each channel)
            ff.close()
            if kline == 0 : 
                if os.path.isfile(filename) :
                    os.remove(filename)
   

def sta_dict_to_pycorr_per_year_and_network(sta_dict,in_fdsn) :
    ''' here we sort the stations we found per year and network, i.e we create one file per year/network
        if a station has two datacenter, we output it twice'''

    out_dir = 'stations_per_year_network'
    if not os.path.isdir(out_dir) :
        os.mkdir(out_dir)

    year1 = in_fdsn['start'].split('-')[0]
    year2 = in_fdsn['end'].split('-')[0]

    for iyear in np.arange(int(year1),int(year2)+1) : # loop on year 
        for inet in sta_dict :                                 # loop on the network/dc/station/ch 
            filename = out_dir+'/'+str(iyear)+'_'+inet+'.txt'
            ff = open(filename, "w") # create the output file 
            kline=0 
            for idc in sta_dict[inet] : 
                for ista in sta_dict[inet][idc] :
                    for ich in sta_dict[inet][idc][ista] :
                        csta = sta_dict[inet][idc][ista][ich]
                        if int(csta['start'].split('-')[0]) > iyear : # si la stationa demarre apres l'anne courrante
                            continue 
                        if int(csta['end'].split('-')[0]) < iyear : # si la station s'est arrete avant l'anne courrante
                            continue 

                        cline = csta['dc']+'   '+csta['net']+'   '+csta['sta'] +'   '+csta['loc']+'   ' + csta['lat']
                        cline = cline + '   '+csta['lon']+'   '+csta['elev']+'   '+csta['depth']+'\n'
                        ff.write(cline)
                        kline = kline+1
                        break # we output only the station one time (and not for each channel)
            ff.close()
            if kline == 0 : 
                if os.path.isfile(filename) :
                    os.remove(filename)



#def sta_dict_to_pycorr_per_year_and_network_replacing_IRIS_ORFEUS(sta_dict,in_fdsn) :
#   ''' here we sort the stations we found per year and network, i.e we create one file per year/network
#        if a station has two datacenter, we output it twice'''
#
#   out_dir = 'stations_per_year_network_replacing_IRIS_ORFEUS_NOA'
#  if not os.path.isdir(out_dir) :
#        os.mkdir(out_dir)
#
#   year1 = in_fdsn['start'].split('-')[0]
#  year2 = in_fdsn['end'].split('-')[0]
#
#    for iyear in np.arange(int(year1),int(year2)+1) : # loop on year 
#        for inet in sta_dict :   # loop on the network/dc/station/ch 
#            filename  = out_dir+'/'+str(iyear)+'_'+inet+'.txt'
#            ff = open(filename, "w") # create the output file 
#            kline = 0 
#            dc_list = list(sta_dict[inet].keys())
#            for idc in sta_dict[inet] : 
#                for ista in sta_dict[inet][idc] :
#                    for ich in sta_dict[inet][idc][ista] :
#                        csta = sta_dict[inet][idc][ista][ich]
#                        if int(csta['start'].split('-')[0]) > iyear : # si la stationa demarre apres l'anne courrante
#                            continue 
#                        if int(csta['end'].split('-')[0]) < iyear : # si la station s'est arrete avant l'anne courrante
#                            continue 
#                        if idc == 'IRIS' : #special case for IRIS !!
#                            for idc2 in ['ORFEUS','NOA'] : 
#                                already_here = 0
#                                if idc2 in sta_dict[inet] :
#                                    if ista in sta_dict[inet][idc2] :
#                                        for ich2 in sta_dict[inet][idc2][ista] :
#                                            if int(sta_dict[inet][idc2][ista][ich2]['start'].split('-')[0]) <= iyear #: 
#                                                if int(sta_dict[inet][idc2][ista][ich2]['end'].split('-')[0]) >= #iyear : 
#                                                    already_here=1 
#                                if already_here ==0 :                                         
#                                    cline = idc2+'   '+csta['net']+'   '+csta['sta'] +'   '+csta['loc']+'   ' + #csta['lat']
#                                    cline = cline + '   '+csta['lon']+'   '+csta['elev']+'   '+csta['depth']+'\n'
#                                    ff.write(cline)
#                                    kline=kline+1
#
#                        cline = csta['dc']+'   '+csta['net']+'   '+csta['sta'] +'   '+csta['loc']+'   ' + csta['lat']
#                        cline = cline + '   '+csta['lon']+'   '+csta['elev']+'   '+csta['depth']+'\n'
#                        ff.write(cline)
#                        kline = kline+1
#                        break # we output only the station one time (and not for each channel)
#            ff.close()
#            if kline == 0 : 
#                if os.path.isfile(filename) :
#                    os.remove(filename)
#
   

def sta_dict_to_pycorr_per_year(sta_dict,in_fdsn) :
    ''' here we sort the stations we found per year, i.e we create one file per year
        if a station has two datacenter, we output it twice'''
    out_dir = 'stations_per_year'
    if not os.path.isdir(out_dir) :
        os.mkdir(out_dir)

    year1 = in_fdsn['start'].split('-')[0]
    year2 = in_fdsn['end'].split('-')[0]
   
    for iyear in np.arange(int(year1),int(year2)+1) : # loop on year 
        filename = out_dir+'/'+str(iyear)+'_stations.txt'
        ff = open(filename, "w") # create the output file 
        kline = 0 
        for inet in sta_dict :                                 # loop on the network/dc/station/ch 
            for idc in sta_dict[inet] : 
                for ista in sta_dict[inet][idc] :
                    is_bb = False 
                    is_sp = False
                    is_acc= False 
                    

                    for ich in sta_dict[inet][idc][ista] :
                        csta = sta_dict[inet][idc][ista][ich]
                        if int(csta['start'].split('-')[0]) > iyear : # si la stationa demarre apres l'anne courrante
                            continue
                        #try :
                        #    int(csta['end'].split('-')[0])
                        #    print(csta['end'])
                        #except : 
                        #    print(csta['end'])
                        #    ipdb.set_trace() 
                        if int(csta['end'].split('-')[0]) < iyear : # si la station s'est arrete avant l'anne courrante
                            continue 

                        cline = csta['dc']+'   '+csta['net']+'   '+csta['sta'] +'   '+csta['loc']+'   ' + csta['lat']
                        cline = cline + '   '+csta['lon']+'   '+csta['elev']+'   '+csta['depth']+'\n'
                        ff.write(cline)
                        break # we output only the station one time (and not for each channel)

    ff.close()
    if kline == 0 : 
        if os.path.isfile(filename) :
            os.remove(filename)



def sta_info_to_dict(sta_info, in_) :
    sta={}

    for iline in sta_info :
        if iline[0:10]== '#Network |' : 
            continue 
        elif len(iline)==0 :
          continue 
        elif iline[0:11]=='#DATACENTER' :
          dc = iline.split(',')[0].split('=')[-1]
          dc = clean_dcname(dc) 
        else :
            csta=iline.split('|')
            lat = float(csta[4])
            lon = float(csta[5])
            if lat < in_['lat'][0] : continue 
            if lat > in_['lat'][1] : continue 
            if lon < in_['lon'][0] : continue 
            if lon > in_['lon'][1] : continue 

            net = csta[0]
            if net not in sta :
                sta[net]={}
            if dc not in sta[net] :
                sta[net][dc]={}

            if len(csta[2])==0 : # if the location code is not set
                csta[2]='--'     # put this ;) 
            sta_id = csta[0]+'_'+csta[1]+'_'+csta[2] # FR.ISO.00
            if sta_id not in sta[net][dc] :
                sta[net][dc][sta_id]={}

            ch = csta[3] 
            if ch not in sta[net][dc][sta_id] : 
                sta[net][dc][sta_id][ch]={}


            sta[net][dc][sta_id][ch]['dc']=dc 
            sta[net][dc][sta_id][ch]['net']=csta[0]
            sta[net][dc][sta_id][ch]['sta']=csta[1]
            sta[net][dc][sta_id][ch]['loc']=csta[2]
            sta[net][dc][sta_id][ch]['ch'] =csta[3]
            sta[net][dc][sta_id][ch]['lat']=csta[4]
            sta[net][dc][sta_id][ch]['lon']=csta[5]
            sta[net][dc][sta_id][ch]['elev']=csta[6]
            sta[net][dc][sta_id][ch]['depth']=csta[7]
            sta[net][dc][sta_id][ch]['azimuth']=csta[8]
            sta[net][dc][sta_id][ch]['dip']=csta[9]
            sta[net][dc][sta_id][ch]['sensor']=csta[10]
            sta[net][dc][sta_id][ch]['start']=csta[15]
            sta[net][dc][sta_id][ch]['end']=csta[16]

            # complete the information if some are missing (extremely rare)
            csta2 = sta[net][dc][sta_id][ch] 
            for ifield in ['depth','elev','azimuth'] :
                if len(csta2[ifield]) == 0 : 
                    csta2[ifield]='0.0'
            for ifield in ['sensor'] : #,'start','end'] :
                if len(csta2[ifield]) == 0 : csta2[ifield]='unknown'
            if len(csta2['dip']) == 0 : 
                csta2['disp']='-999' 
    return sta

def sta_dict_filter_channel(sta_dict,in_) :
#''' if in_['sp_but_no_bb']== True : 
#        we keep only stations that have 1- a short period channel [EH,SH,MH] and 2) no broadband [BH,HH]
#        broadband stations and accelerometers are discarded
#'''

    # if the user did not ask any filter we return : 
    if in_['sp_but_no_bb'] == False and in_['acc_but_no_vel']==False :
        return sta_dict 
    dd.dispc(' ','b','n')
    dd.dispc(  'filtering channel as requested by user','c','b')

    if in_['sp_but_no_bb']== True : 
        dd.dispc('  in_[sp_but_no_bb] = True, so we keep only stations if :','c','n') 
        dd.dispc('    1) there is a short period channel [EH,SH,MH]','c','n')
        dd.dispc('    2) there is no broadband channel [BH,HH]','c','n')
        dd.dispc('  in any case broadband and accelerometers channel are discarded','c','n')
        dd.dispc('  this options makes only sense if in_[channel] is set to *Z','c','n')

    if in_['acc_but_no_vel']== True : 
        dd.dispc('  in_[acc_but_no_vel] = True, so we keep only stations if :','c','n') 
        dd.dispc('    1) there is an accelerometer channel [*N*]','c','n')
        dd.dispc('    2) there is no velocimetric channel [*H*,*H*]','c','n')
        dd.dispc('  in any case SP and BB channels are discarded','c','n')
        dd.dispc('  this options makes only sense if in_[channel] is set to * or *Z','c','n')
    dd.dispc('','b','n')


    # scann each station and to list all channels belonging to any location code 
    for inet in sta_dict : 
        for idc in sta_dict[inet] : 
            # we are in sta['FR']['RESIF'] : we build a dictionnary of channels
            per_sta = {}
            #per_sta['id_']   = []
            #per_sta['ch']    = [] 
            cnet = sta_dict[inet][idc]
            for ista in cnet : 
                id_short = ista.split('_')[0]+'_'+ista.split('_')[1]
                if not id_short in per_sta :
                    per_sta[id_short]={}
                    per_sta[id_short]['id_']=[]
                    per_sta[id_short]['ch'] =[]
                for ich in cnet[ista] :
                   per_sta[id_short]['id_'].append(ista)
                   per_sta[id_short]['ch'].append(ich)

            # now determine which channel is an BB, SP or ACC
            for ista in per_sta :
                csta = per_sta[ista]
                nch = len(csta['ch'])
                csta['is_bb'] = np.zeros(nch)
                csta['is_sp'] = np.zeros(nch)
                csta['is_acc'] = np.zeros(nch)
                k=-1 
                for ich in csta['ch'] : 
                    k=k+1
                    if ich[0:2] == 'HH' or ich[0:2] =='BH' or ich[0:2] =='LH' or ich[0:2] =='VH' :
                        csta['is_bb'][k]= True
                    elif ich[0:2] == 'MH' or ich[0:2] == 'SH' or ich[0:2] == 'EH' :
                        csta['is_sp'][k]=True
                    elif ich[1:2] == 'N' :
                        csta['is_acc'][k]=True


                if in_['sp_but_no_bb'] == True :    
                    if csta['is_bb'].sum() > 0 :           # we have some BB so we have to remove all SP
                        I=np.where(csta['is_sp']==True)[0] # get indice of all SP
                        for iI in I :                      # remove them from sta_dict 
                            sta_dict[inet][idc][csta['id_'][iI]].pop(csta['ch'][iI])
                    #in any case remove broadband and accelerometers : 
                    I=np.where(csta['is_bb']+csta['is_acc']>0)[0] # get indice of all BB or ACC
                    for iI in I :               
                        sta_dict[inet][idc][csta['id_'][iI]].pop(csta['ch'][iI]) # remove them from dict

                if in_['acc_but_no_vel'] == True :    
                    if csta['is_bb'].sum() + csta['is_sp'].sum() > 0 : # we have some BB OR SP so we have to rm all ACC
                        I=np.where(csta['is_acc']==True)[0] # get indice of all SP
                        for iI in I :                      # remove them from sta_dict 
                            sta_dict[inet][idc][csta['id_'][iI]].pop(csta['ch'][iI])
                    #in any case remove velocimeters :
                    I=np.where(csta['is_bb']+csta['is_sp']>0)[0] # get indice of all BB or ACC
                    for iI in I :               
                        sta_dict[inet][idc][csta['id_'][iI]].pop(csta['ch'][iI]) # remove them from dict


    return sta_dict 


def sta_dict_remove_empty_net_dc_sta(sta_dict) : 
    # first removing empty stations :
    for inet in sta_dict : 
        for idc in sta_dict[inet] : 
            sta_list = []
            for ista in sta_dict[inet][idc] :
                if not sta_dict[inet][idc][ista] : 
                    sta_list.append(ista)
            for ista in set(sta_list) :
                sta_dict[inet][idc].pop(ista)

    # then removing empty dc :
    for inet in sta_dict : 
        dc_list =[]
        for idc in sta_dict[inet] : 
            if not sta_dict[inet][idc]: 
                    dc_list.append(idc)
        for idc in set(dc_list) :
            sta_dict[inet].pop(idc)

    #finally remove empty network: 
    net_list =[]
    for inet in sta_dict : 
        if not sta_dict[inet]: 
            net_list.append(inet)
    for inet in set(net_list) :
        sta_dict.pop(inet)
    return sta_dict

def sta_dict_count(sta_dict,str_) : 
    net = [] 
    dc  = [] 
    id_ = [] 
    ch  = [] 
    for inet in sta_dict : 
        net.append(inet)
        for idc in sta_dict[inet] : 
            dc.append(idc)
            for ista in sta_dict[inet][idc] :
                id_.append(inet+'_'+'ista')
                for ich in sta_dict[inet][idc][ista] :
                    ch.append(ich) 

    dd.dispc(str_,'c','b')
    dd.dispc('    '+str(len(net))+' networks','c','n')
    dd.dispc('    '+str(len(id_))+ ' stations','c','n')
    dd.dispc('    '+str(len(set(ch))) + ' different kind of channels','c','n')
    dd.dispc('    '+str(len(dc))     + ' data centers','c','n')

def clean_dcname(dc) :
    if dc=='IRISDMC' :dc = 'IRIS'
    if dc=='SED'     :dc = 'ETH'
    if dc=='GEOFON'  :dc = 'GFZ'
    if dc=='USPSC'   :dc = 'USP'
    return dc 

def print_documentation(in_iris,in_) : 

    dd.dispc(' ','w','n')
    dd.dispc('to get all velocimeters use :','w','b')
    dd.dispc("in_iris['channel']         = '*HZ'",'w','n')
    dd.dispc("in_['sp_but_no_bb']        = False ",'w','n')
    dd.dispc("in_['acc_but_no_vel']      = False ",'w','n')
    dd.dispc("get_sta.find_stations(in_iris,in_)",'w','n')
    dd.dispc('','w','b')
    dd.dispc('to get all broadband velocimeters :','w','b')
    dd.dispc("in_iris['channel']         = 'HHZ,BHZ'",'w','n')
    dd.dispc("in_['sp_but_no_bb']        = False ",'w','n')
    dd.dispc("in_['acc_but_no_vel']      = False",'w','n')
    dd.dispc("get_sta.find_stations(in_iris,in_)",'w','n')
    dd.dispc('','w','b')
    dd.dispc('to get short periods velocimeters for stations having no broadband channel :','w','b')
    dd.dispc("in_iris['channel']         = '*HZ'",'w','n')
    dd.dispc("in_['sp_but_no_bb']        = True ",'w','n')
    dd.dispc("in_['acc_but_no_vel']      = False ",'w','n')
    dd.dispc("get_sta.find_stations(in_iris,in_)",'w','n')
    dd.dispc('','w','n')
    dd.dispc('to get accelerometers only for stations having no velocimeter :','w','b')
    dd.dispc("in_iris['channel']         = '*Z'",'w','n')
    dd.dispc("in_['sp_but_no_bb']        = False",'w','n')
    dd.dispc("in_['acc_but_no_vel']      = True",'w','n')
    dd.dispc("get_sta.find_stations(in_iris,in_)",'w','n')
    dd.dispc('','w','n')



##########################################################################################
##########################################################################################
#########           #####         ####               #####           #####################
#########           ####          ####               ####     ############################
#########     ##########     ##############     #########     ############################
#########     #########     ###############     ##########        ########################
#########        ######     ###############     ##############       #####################
#########     #########      ##############     ##################    ####################
#########     ##########          #########     #################     ####################
#########     ###########         #########     ##########           #####################
##########################################################################################
##########################################################################################



#---------------------------------------------------------------------------------------
#-------------------- STAIONS SUBFUNCTIONSFUNCTIONS ----------------------------------
#---------------------------------------------------------------------------------------


def mmap(lon=[2.39,5.72],lat=[47.08,45.18],filename='stations.png',info='station map'):
    stamen_terrain = Stamen('terrain-background')
    ax = plt.axes(projection=stamen_terrain.crs)
    minlat = -80
    maxlat = 80
    minlon = -180
    maxlon = 180
    dl     = 0.1
    r      = int(locations2degrees(min(lat),min(lon),max(lat),max(lon))*111)
    if r < 10000:
        if min(lat) - dl > minlat: minlat = min(lat) - dl
        if max(lat) + dl < maxlat: maxlat = max(lat) + dl
        if min(lon) - dl > minlon: minlon = min(lon) - dl
        if max(lon) + dl < maxlon: maxlon = max(lon) + dl
    ax.set_extent([minlon, maxlon, minlat, maxlat])
    r = int(locations2degrees(minlat,minlon,maxlat,maxlon)*111)
    if r <= 1000 :
        res = 8 
        ax.coastlines('10m')
    if r > 1000 and r <= 5000:
        res = 6
        ax.coastlines('50m')    
    if r > 5000 and r <= 10000: 
        res = 4
        ax.coastlines('110m')
    if r > 10000: 
        res = 2
        ax.coastlines('110m')
    ax.add_image(stamen_terrain,res)
    gl = ax.gridlines(draw_labels=True,linewidth=0.5, color='gray', alpha=0.5)
    gl.xlabels_top   = False
    gl.ylabels_right = False    
    plt.scatter(lon, lat, c='red', s=20,transform=ccrs.Geodetic(),cmap='plasma_r',alpha=0.5)
    plt.savefig(filename,dpi=150, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches='tight', pad_inches=0.1,
        frameon=None)
    #plt.show()
    plt.close()







