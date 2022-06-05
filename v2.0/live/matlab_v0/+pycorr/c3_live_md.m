%---------------------------------------------------------------------
% Typical usage : 
% h=pycorr.c3_live_md('C3__05_c3__05-10s__coda_10dp-1000s_normenv')
% h.p_map()
%
%--------------------------------------------------------------------
% PLOTTING METHODS 
%   p_station        : plot all virtual source & sta + sta in a != color
%     -h.p_map 
%   p_station_pair   : plot all virtual source & sta + sta pair in a != color
%     -h.p_map
%   p_map            : plot  all virtual source + station 
%
%----------------------------------------------------------------------
% UTILITY METHODS USED BY C3 LIVE : 
%   get_path_info : for a given path returns a structure with the h5 file path, h5 filename, 
%                   a suffix + prefix useful to determine the path plot title
%                   these information are used by the C1_trace constructor. 
%  
%      - get_prefix_from_dtype    : return '/c3','/c1','c3_daily_stacked' 
%      - get_filename_from_indice :
%
% ----------------------------------------------------------------------------
% UTILITY METHODS : 
%   get_station_list : get the list of station 
%
%------------------------------------------------------------------------
% live_c3 : generate a db.mat file and read it 
%    => build_db_file 
%       1.read metadata of each C3*/xcorr_set_0000.h5 file 
%       2.rearange and cct them 
%       3.compute dist,az,baz of each station pair 
%       4.remember in which file each path is stored (h.file_I)
%      => read_metadata(h5_file) : 
%        return metadata of a single h5 file 
%      => celltrim : used to trim h.id 
%-----------------------------------------------------------------------
classdef c3_live_md < handle 
	properties 
		id    = [];
		lon   = [];
		lat   = [];
		elev  = [];
		dist  = [];
		az    = []; 
		baz   = [];
		t     = []; 
		tau   = []; 
		t_c1a = []; 
		t_c1b = [];
		date1 = []; 
		date2 = [];
		file  = {}; 
		file_I= []; 
		in  = struct();
		md_c= struct();
		vs  = struct();
	end
	%---------------------------------------
	%  plotting metadata 
	%---------------------------------------
	methods 
		function p_station(h,ksta,varargin)
	  	in.color='r'     ;  
		  in.mark='^'      ; 
  		in.msize=10      ;
  		in.border='k'    ;
  		in=lang.parse_options(in,varargin); 
	  	sta=h.get_station_list;
	  	h.p_map(struct(),struct(),struct('color',[0.5 0.5 0.5],'mark','o'))
	  	map.station(sta.lon(ksta),sta.lat(ksta),in,'mark'); 
	  end

		function p_station_pair(h,kpath)
	  	h.p_map(struct(),struct(),struct('color',[0.5 0.5 0.5],'mark','o','msize',4))
	  	for ipath = kpath 
	  		map.station(h.lon(:,ipath),h.lat(:,ipath),'color','r','mark','^','msize',8)
	  	end
	  end

		function p_map(h,in_map,in_background,in_sta)
	  	if ~exist('in_map') ; in_map = struct() ; end 
	  	if ~exist('in_sta') ; in_sta = struct('mark','^') ; end 
	  	if ~exist('in_background') ; in_background = struct() ; end
	  	sta=h.get_station_list;
	  	I=find(sta.lon~=0 & sta.lat ~=0) ;
	  	sta.lon = sta.lon(I); 
	  	sta.lat = sta.lat(I); 
	  	sta.sta = sta.sta(I);
	  	sta.elev= sta.elev(I);
	  	in.lon(1) = min(sta.lon);%-1; 
	  	in.lon(2) = max(sta.lon);%+1; 
	  	in.lat(1) = min(sta.lat);%-1;
	  	in.lat(2) = max(sta.lat);%+1;
	  	in=lang.parse_options(in,in_map);
	  	map.background(in.lon,in.lat,in_background)
	  	map.station(h.vs.lon,h.vs.lat,'color','w','msize',4,'mark','^');
	  	map.station(sta.lon,sta.lat ,in_sta); %'color','r','msize',4,'mark','o')
	  end

	  function sta=get_station_list(h)
	  	[sta_ Ista]=unique(h.id(:));
      sta.lon=double(h.lon(Ista));
      sta.lat=double(h.lat(Ista));
      sta.elev=double(h.elev(Ista));
      sta.sta=h.id(Ista);
	  end
	end
	%--------------------------------------------------------------------
	% UTILITY METHODS 
	%----------------------------------------------------------------------
	methods 
		function path_info = get_path_info(h,kpath,c3_type,dtype,kcmp) 
			prefix = h.get_prefix_from_dtype(dtype)                  ; 
			if dtype == 0 % c1 
				path_  = ['/c1/',h.id{1,kpath},'/',h.id{2,kpath},'/ZZ'];
				suffix = '' ;
				cmp    ='ZZ';
			else 		
				path_    = [prefix,'/',h.id{1,kpath},'/',h.id{2,kpath}];
				path_    = [path_,'/',h.in.cmp{kcmp}]                  ;
				path_    = [path_, '/',h.in.c3_type{c3_type}]          ;
				suffix   = h.in.c3_type{c3_type}                       ;
				cmp      = h.in.cmp{kcmp}                              ;
			end
			fname = h.get_filename_from_indice(kpath) ;
			path_info.path_ = path_ ; 
			path_info.prefix= prefix; 
			path_info.suffix= suffix; 
			path_info.fname = fname ; 
			path_info.cmp   =cmp    ;
		end	  	

		function prefix=get_prefix_from_dtype(h,dtype)
			if dtype == 0 ; 
				prefix = '/c1' ; 
			elseif dtype == 1 
				prefix = '/c3' ;
			elseif dtype == 2 
				prefix = '/c3_daily_stacked';
			elseif dtype==3 
				prefix ='/c3_mat';
			elseif dtype ==4
				prefix='/c3_daily_mat_stacked';
			elseif dtype==5 
				prefix='/c3_daily';
			end
		end

	  function [filename]=get_filename_from_indice(h,kpath)
		 	% return the filename of the correlation #I 
			I_file=find(h.file_I(2,:) >= kpath,1,'first') ;
	 	  filename=h.file{I_file}                       ;
		end 
	end

	methods 

		%constructor : generate a db.mat file it does not exist and load it
		function h=c3_live_md(in_dir)
			db_file=[in_dir,'/db.mat']; 
			if ~exist(db_file,'file')
				dispc('  creating new db file','y','b')
				db = build_db_file(in_dir)  ;
				save(db_file,'-struct','db');
			end
			dispc('  loading existing db file','y','b');
			h=struct_.copy_field_to(load(db_file),h);
		end
	end
end

%--------------------------------------------------------------
% FUNCTION OUTSIDE THE CLASS USED TO PREPARE THE DB.mat file 
%-------------------------------------------------------------
function db=build_db_file(in_dir)
	%1. read the metadata of each file 
	set_list = os.dir([in_dir,'/xcorr_set*h5']);
	for iset = set_list 
		cset=iset{1}; 
		if ~exist('md_all')
			md_all=read_metadata(cset)        ;
		else 
			md_all(end+1)=read_metadata(cset) ;
		end
	end
	%2. copy them in to db 
	%2.1 get the number of path 
	nfile = numel(md_all);
	db = struct(); 
	npath =0 ;
	for ifile =1 : nfile 
		npath = npath + numel(md_all(ifile).md.id(1,:)); 
	end 
	%2.2 initialize and copy data into db : 
	db.lat = zeros(2,npath); 
	db.lon = zeros(2,npath); 
	db.id  = cell(2,npath) ; 
	db.elev= zeros(2,npath); 
	I1=1;
	for ifile =1 : nfile 
		npath_in_this_file = numel(md_all(ifile).md.id(1,:)); 
		I2 = I1 + npath_in_this_file-1          ;
		db.lat(:,I1:I2) = md_all(ifile).md.lat  ; 
		db.lon(:,I1:I2) = md_all(ifile).md.lon  ; 
		db.id(:,I1:I2) = md_all(ifile).md.id    ;
		db.elev(:,I1:I2) = md_all(ifile).md.elev;
		I1=I2+1;
	end
	%2.3 take care of md_c and in_ , vs
	db.vs   = md_all(1).vs         ;
	db.in   = md_all(1).in         ;
	db.md_c = md_all(1).md_c       ;
	db.t    = md_all(1).md_c.t     ;
	db.t_c1a= md_all(1).md_c.t_c1a ;
	db.t_c1b= md_all(1).md_c.t_c1b ;
	db.tau  = md_all(1).md_c.tau   ; 
	if isfield(md_all(1).md_c,'date1_ordinal')
		db.date1= md_all(1).md_c.date1_ordinal ; 
		db.date2= md_all(1).md_c.date2_ordinal ; 
	else %when not doing monitoring  
		db.date1=datenum(2000,1,1); 
		db.date2=datenum(2000,1,2);
	end
	%2.4 compute distance and azimuth for all couple of station : 
	db.dist = zeros(1,npath);
	db.az   = zeros(1,npath);
	db.baz  = zeros(1,npath);
	for ipath =1 : npath 
		lon1 = db.lon(1,ipath); 
		lon2 = db.lon(2,ipath);
		lat1 = db.lat(1,ipath);
		lat2 = db.lat(2,ipath);
		[db.dist(ipath) db.az(ipath) db.baz(ipath)]=m_idist(lon1,lat1,lon2,lat2);
	end
	db.dist=db.dist/1000; 
	db.dist(find(isnan(db.dist)==1))=0;
	%2.5 need to remember in which file are stored each c3 : 
	db.file  =cell(1,nfile) ;
	db.file_I=zeros(2,nfile);
	I1=1;
	for ifile =1 :nfile 
		db.file(ifile)=set_list(ifile)      ;
		I2=I1+size(md_all(ifile).md.lon,2)-1;
		db.file_I(:,ifile)=[I1 I2]          ;
		I1=I2+1                             ;
	end
	%2.6 convert python dates into matlab 
	%2.7 remove white space at the end of ID name : 
	db.id =cellstr_.trim(db.id);
	%db.id=cellfun(@celltrim,db.id)      ;
end


function toto=read_metadata(h5_file)
	toto = struct()                             ;
	toto.md   = h5.read_group(h5_file,'/md')    ;
	toto.md_c = h5.read_group(h5_file,'/md_c')  ;
	toto.vs   = h5.read_group(h5_file,'/vs')    ; 
	toto.in   = h5.read_group(h5_file,'/in_')   ;
	toto.in.gw= h5.read_group(h5_file,'/in_/gw');
	toto.in.pp= h5.read_group(h5_file,'/in_/pp');
	toto.in.h1= h5.read_group(h5_file,'/in_h1') ;
	toto.in.h2= h5.read_group(h5_file,'/in_h2') ;
end
