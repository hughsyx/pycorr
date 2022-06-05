%this script is useful to check the selection of virtual sources used to convert C2 into C1 : 
%in that case the C3_to_C1 script should be run the "save_mat" option set to true. 

file=os.dir('test/C3_to_C1__dist_15_60km__daz_5__norm__rm_btw_sta/*mat');


for ifile=1:numel(file);
    toto=file{ifile}
    h=load(toto(2:end-1)); 
    if max(h.cc_sum)>0
        clf;
        ha=p.split_map_trace; 
        subplot(ha(1)); 
        lon =[min([h.lon h.vs.lon]) max([h.lon h.vs.lon])];
        lat =[min([h.lat h.vs.lat]) max([h.lat h.vs.lat])];
        map.background(lon,lat,'topo',false);
        map.station(h.vs.lon,h.vs.lat,'color','w')
        map.station(h.vs.lon(h.vs.I_sel+1),h.vs.lat(h.vs.I_sel+1),'color','y')
        map.station(h.lon,h.lat,'color','r')
        %map.station(h.vs.lon(h.vs.I_az),h.vs.lat(h.vs.I_az),'color','c')
        subplot(ha(2))
        hold on;
        plot(h.time,h.cc_pp_I1,'r')
        plot(h.time,h.cc_nn_I1,'b')
        plot(h.time,h.cc_pp_I2,'c')
        plot(h.time,h.cc_nn_I2,'g')
        plot(h.time,s2d.filter(h.cc_sum,1,5,0.04),'k')
        h.vs_info.dist1(h.vs.I_sel+1)
        h.vs_info.dist2(h.vs.I_sel+1)
    end
end
