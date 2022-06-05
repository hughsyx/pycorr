% reading a single reference correlation directly from matlab
% without using the pycorr library : 


kpath = 10  ; % we want to read the trace of the station pair #kpath 
kcmp = 1    ; % for the 1st component 
file = 'C1/xcorr_all/xcorr_CH_CH_0000.h5';

% first we read all correlations metadata of this file : 
h = struct();
h.md   = h5_read_group(file,'/md')  ;
h.md_c = h5_read_group(file,'/md_c');
h.md.id=cellfun(@celltrim,h.md.id)  ; % remove extrat white space 

% now we get the path of this trace kpath within the h5 file : 
path_=['/ref/',h.md.id{1,kpath},'/',h.md.id{2,kpath},'/',h.md_c.cmp{kcmp}];

% reading the trace : 
cc=h5read(file,[path_]);  
plot(h.md_c.t,cc); 
xlim([-100 100]);
title([h.md.id{1,kpath},'-',h.md.id{2,kpath},'  ', h.md_c.cmp{kcmp}]);


function group_data=h5_read_group(h5_filename,group)
    group_data=struct();
    group_info=h5info(h5_filename,group);
    for idset = [1:numel(group_info.Datasets)];
        dname=group_info.Datasets(idset).Name;
	group_data.(dname)=h5read(h5_filename,[group,'/',dname]);
    end	
end

function str_=celltrim(str_)
    str_={str_(find(str_~=char(0)))};
end
