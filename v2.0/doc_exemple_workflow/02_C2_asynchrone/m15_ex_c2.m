% single source 
h=pycorr.live('C3_to_C1__single_source__dist_0_1000km__daz_10__norm__rm_btw_sta');

clf; hold on;

npath = numel(h.dist);
npr=0;
for ipath=1 : npath 
    cc=h.read_cc(ipath,1,5);
    if max(cc.cc) >0
        clf; 
        h.p_cc(ipath,1,5)
        p.png('large4','dpi',100,'tag',['_#',num2str(ipath)]);
        npr= npr+1; 
    end
    if npr > 10;
        break
    end
end

return
