
% these codes requires the m_map packages to be installed 

% add the matlab library manually but it is of course better to add
% these path in your. bashrc file :
addpath('matlab')
addpath('matlab/library/functions')
addpath('matlab/library/packages_1')

h=pycorr.live('C1/xcorr_all'); 


% plot station map. Install m_map etopo file if you want it to look
% great 
h.p_map 



% retrieve the list of station :
h.get_station_list

% plot the 1st path filtered btw 5-10s : 
clf ; h.p_cc(1,5,10);

% plot the 1st path filtered btw 5-10s, 1st component (ZZ) : 
clf ; h.p_cc(1,5,10,'cmp',1);


% plot the path 1,2,3  :
clf ; h.p_cc([1 2 3],5,10)

% plot the path 2 filtered in the 5-10s and 10-20s : 
clf ; h.p_cc(2,[5 10],[10 20])

% read the 1st pathm in the 5-10s period band, 1st component (ZZ) 
cc = h.read_single_cc(1,5,10,1);

% read all CC computed btw station 1 and the rest of the network:
% cc=h.read_cc_sta_vs_others(1,5,10)  

