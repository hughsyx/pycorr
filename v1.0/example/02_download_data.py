import m_pycorr.m10_get_data as get_data

############### 
# Standard
############### 
in_ = {}
in_['path']='./data_1.0hz/daily/glob/2018'
in_['mode']= 0

############### RUN
get_data.download(in_)


