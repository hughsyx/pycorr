% single source 
h=pycorr.live('C3_to_C1__dist_5_60km__daz_5__norm__rm_btw_sta')

clf; hold on;
h.p_cc_vs_dist(1,5,'keep_trace_with_data',1,'npath',inf)

xlim([-50 50])
ylim([2 16]);
line([0  100],[0  100*1],'color','r')
line([0  100],[0  100*0.7],'color','b')
legend(' ','1 km/s','0.7 km/s')


aa=strsplit(h.file{1},'/');
aa=aa{1};
net1=strsplit(h.id{1,1},'.'); 
net1=net1{1};
%..
net2=strsplit(h.id{2,1},'.'); 
net2=net2{1};

title({['C2 profile 1',num2str(net1),'- profile 2',num2str(net2)],aa});
p.png('big','tag','all_dist');

ylim([2 14]);
p.png('big','tag','less_20');
