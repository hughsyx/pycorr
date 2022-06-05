% dvv_med_all : median dt/t for one path, all coda window 
% dvv_med_all_std : the corresponding std 
% dvv_med_sel : median dt/t on selected coda window 
% dvv_med_std : corresponding standard deviation 
% dvv_med_less_than_one_std : 


% mapping entre les dvv et les cc : 
% h.dist(h.dvv.I) <=> h.dvv.dist 
% h.dist(h.dvv.I(N)) = h.dvv(N)
classdef dvv_doublet < pycorr.live 
  properties 
    dvv = struct() 
    
  end
  methods 

  end
  methods
    function  dvv_plot(h,I_path,varargin); 
      % loop on date to remove path with no data : 
      npath = numel(I_path)
      ndate = numel(h.dvv.date); 
      dvv = zeros(1,ndate); 
      for idate = 1 : ndate 
        dvv_date = h.dvv.dvv_lg(idate,I_path)   ; 
        I = find(dvv_date ~=0)       ; 
        dvv(idate) = median(dvv_date(I));
      end
      plot(h.dvv.date,dvv)
      datetick('x')
    end
    function I = dvv_select_path(h,varargin)
      %  select station pair by on, lat, and if the dv/v measured each day has a low standard deviation :
      %    for each path, we computed the std of the dv/v measured on different coda window 
      %    we keep only measurements the median(standard deviation vs day) < in.std. 
      %    this ensure that we use the same path every day. 
      in.lon = [-inf inf] ; % longitude window
      in.lat = [-inf inf] ; % latitude window 
      in.std = inf; % standard deviation threshold 
      in = lang.parse_options(in,varargin);
      
      I = [] ; 
      % get the coordinate of the path in the correct order : 
      lon = h.lon(:,h.dvv.I) ; 
      lat = h.lat(:,h.dvv.I) ;
      % look for station in the right area : 
      I1 = find(lon(1,:) <= in.lon(2) & lon(1,:) >= in.lon(1));
      I2 = find(lon(2,:) <= in.lon(2) & lon(2,:) >= in.lon(1));
      I3 = find(lat(1,:) <= in.lat(2) & lat(1,:) >= in.lat(1));
      I4 = find(lat(2,:) <= in.lat(2) & lat(2,:) >= in.lat(1));
      I = intersect(I1,I2) ; 
      I = intersect(I,I3); 
      I = intersect(I,I4); 
      I = unique(I); 
      % now get accurate dvv/v measurements we want to get the same path for all days: 
      I_std = find(median(h.dvv.dvv_med_all_std) < in.std) ; 
      I = intersect(I,I_std); 
    end
  end 
  methods
    function h=dvv_doublet(pydb_file)
      h = h@pycorr.live(pydb_file)                                        ;
      h.dvv=h5.read_group(pydb_file,'/dvv_doublet')                       ;
      %h.dvv.dist  = h5.tree_to_matrix(pydb_file,'/dvv_doublet/dist')      ;
      [~, dset_path]= h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_lg',false) ;
      h.dvv.dvv_lg  = h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_lg')  ;
      h.dvv.dvv_lgb  = h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_lgb');
      h.dvv.dvv_lgr  = h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_lgr');
      h.dvv.dvv_lgp  = h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_lgp');
      h.dvv.dvv_lgstd= h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_lgstd');
      h.dvv.dvv_med_all    = h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_med_all')     ;
      h.dvv.dvv_med_all_std= h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_med_all_std') ;
      h.dvv.dvv_med_sel    = h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_med_sel')     ;
      h.dvv.dvv_med_std= h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_med_std') ;
      h.dvv.dvv_med_sel_less_one_std= h5.tree_to_matrix(pydb_file,'/dvv_doublet/dvv_med_sel_less_one_std') ;
      h.dvv.in          = h5.read_group(pydb_file,'/dvv_doublet/in_')  ;
      h.dvv.in_ml       = h5.read_group(pydb_file,'/dvv_doublet/in_ml');
      h.dvv.md_c        = h5.read_group(pydb_file,'/dvv_doublet/md_c') ;
      h.dvv.date1 = h.dvv.md_c.date1+366 ;
      h.dvv.date2 = h.dvv.md_c.date2+366 ;
      h.dvv.date  = h.dvv.md_c.date +366 ;
      % Establish the mapping between correlation metadata and computed dispersion curve :
      % useful to display the correlation corresponding to a disp curve :
      h_id=cell_.merge(h.id(1,:),h.id(2,:),'/');
      npath = numel(dset_path);
      h.dvv.I =zeros(1,npath);
      h.dvv.id = dset_path;
      for ipath = 1 : numel(dset_path)
        try 
          h.dvv.I(ipath)=strmatch(dset_path{ipath},h_id);
        catch 
          disp('xx')
          keyboard
        end
      end 
    end
  end
end