
classdef live_md < handle 
    properties
        ref_nstack=[] ; % number of days correlated [path,]
    	date1=[] ; % datenum of the 1st point correlated for each days : 
        date2=[] ; % datenum of the last point correlated for each day
    	lon = [] ; % longitude of each pair of station
    	lat = [] ; % latitude  of each pair of station
        elev= [] ; % elevation of each pair of station 
        depth=[] ; 
    	id={}    ; % id os each station pair  la station
    	dist= [] ; % distance between each pair of station 
        az  = [] ; % azimuth of each pair of station 
        baz = [] ; % baz of each pair of station 
        in  = struct(); % user input 
    	version =[];    % version number of the hdf5 file format
    	cmp={}    ;
    	t=[]      ;     % correlation time vector [s]
    	tau=1     ;     % correlation sample interval [s]
    	file={}   ;     % cc.file=relative path file where each station-pair  are stored.
    	file_I=[] ;
    end
    methods 
        function sisp_export_station_list(h)
            sta = h.get_station_list()
            fname = 'STACoord.dat'; 
            ff=fopen(fname,'w');
            for ista = 1 : numel(sta.lon);
                fprintf(ff,'%s 00 %3.8f %3.8f 0\n',sta.sta{ista},sta.lat(ista),sta.lon(ista)) ; 
            end
            fclose(ff);
        end
        %-----------------------------------------------
        %        plotting functions 
        %-----------------------------------------------
        function ha=p_availability(h,varargin)
            in.gmap = true ;
            in.pmap = true ;
            in.pbar = true ;
            in = lang.parse_options(in,varargin);
            nplot= 0 ; 
            if in.pmap 
                nplot = nplot+ 1 ; 
            end
            if in.pbar 
                nplot = nplot + 1; 
            end 
            %assume all autocorr have been computed. we get the number of day from the nday used 
            %to compute the autcorr (ref_nstack). 
            %get autocorr :
            I_ac=find(h.dist==0) ;
            %get station list :
            sta=h.get_station_list; 
            nsta=numel(sta.lon)   ;
            sta.nday=zeros(1,nsta);
            %look at indice with the station list of all autocorr and the number of day:
            for ista = I_ac 
                I_sta=strmatch(h.id{1,ista},sta.sta);
                sta.nday(I_sta)=h.ref_nstack(1,ista);
            end
            %plot map with colored stations : 
            kplot=1;
            if in.pmap 
                ha(1)=subplot(nplot,1,kplot)
                I=find(sta.lon~=0 | sta.lat~=0);
                if ~in.gmap 
                    h.p_map(struct(),struct('topo',false),struct('mark','^'));	  	
                    map.station_color(sta.lon(I),sta.lat(I),sta.nday,struct('mark','^'));			
                else 
                    gmap.station_color(sta.lon(I),sta.lat(I),sta.nday,struct('mark','^'));
                end
                kplot = kplot + 1 ;
                title('number of days recorded','fontSize',14,'fontweight','bold')
            end
            %plot results as bar : 
            if in.pbar 
                ha(2)=subplot(nplot,1,kplot)
                bar(sta.nday);
                ylabel('nday');
                set(gca,'Xtick',[1:nsta])    ;
                set(gca,'XtickLabel',sta.sta);
                p.rotate_tick_label(gca,90)  ;
                xlim([0 nsta+1]); 
                titre=['number of days used to compute the autocorr'];
                title(titre,'fontweight','bold','fontSize',14)
            end
        end
        
        function p_station(h,ksta,varargin)
            in.color='r'     ;  %cc
            in.mark='square' ; 
            in.msize=10      ;
            in.border='k'    ;
            in=lang.parse_options(in,varargin); 
            sta=h.get_station_list;
            map.station(sta.lon(ksta),sta.lat(ksta),in); 
        end

        function titre = p_station_pair_gmap(h,kpath,varargin)
            in.col='r'; 
            in.mark='^'; 
            in.msize=12;
            in.type='terrain';
            in.scale =2 ;
            in.sl=0 ;
            in.resize =1 ;
            in.alpha =1;
            in = lang.parse_options(in,varargin);
            col=[0.5  0.5  0.5];
            h.p_gmap('color',col,'msize',in.msize,'type',in.type,'scale',in.scale,'resize',in.resize,'alpha',in.alpha);
            hold on; 
            for ipath=kpath
                plot(h.lon(1,ipath),h.lat(1,ipath),[in.mark,'k'],'MarkerSize',in.msize,'MarkerFaceColor',in.col)
                plot(h.lon(2,ipath),h.lat(2,ipath),[in.mark,'k'],'MarkerSize',in.msize,'MarkerFaceColor',in.col)
            end
            if numel(ipath)==1;
                titre = h.get_title(kpath); 
                title(titre,'fontsize',14,'fontweight','bold');
            end
        end
        function p_station_pair(h,kpath,msize)
            if ~exist('msize') 
                msize=12;
            end
            h.p_map(struct(),struct(),struct('color',[0.5 0.5 0.5],'msize',msize))
            hold on;
            for ipath = kpath 
                map.station(h.lon(:,ipath),h.lat(:,ipath),'color','r','msize',msize)
            end
        end

        function p_gmap(h,varargin)
            in.color='r';
            in.msize=10 ;
            in.type='terrain';
            in.scale =1 ;
            in.sl=0 ;
            in.resize =1 ;
            in.alpha =1;
            in = lang.parse_options(in,varargin);
            sta= h.get_station_list             ;
            I=find(sta.lon~=0 | sta.lat~=0)     ;
            gmap.map(sta.lon(I),sta.lat(I),in)  ;
        end

        function p_map(h,in_map,in_background,in_sta)
            if ~exist('in_map') ; in_map = struct() ; end 
            if ~exist('in_sta') ; in_sta = struct() ; end 
            if ~exist('in_background') ; in_background = struct() ; end
            sta=h.get_station_list;
            I=find(sta.lon~=0 | sta.lat~=0);
            in.lon(1) = min(sta.lon(I));%-1; 
            in.lon(2) = max(sta.lon(I));%+1; 
            in.lat(1) = min(sta.lat(I));%-1;
            in.lat(2) = max(sta.lat(I));%+1;
            in=lang.parse_options(in,in_map);
            map.background(in.lon,in.lat,in_background);
            map.station(sta.lon(I),sta.lat(I) ,in_sta);
        end

        %----------------------------------------------
        %        get correlation indice               
        %----------------------------------------------
        function [I_id_net I_net_id sta_info]=get_indice_path_with_station(h,ksta)
            sta=h.get_station_list; 
            id_ref=sta.sta{ksta}  ;
            %get all correlation id_ref -> network : 
            I_id_net=find(strcmp(h.id(1,:),id_ref)==1); 
            %get all correlation network -> id_ref : 
            I_net_id    =find(strcmp(h.id(2,:),id_ref)==1); 
            %   we need to remove the autocorr since we would get it twice otherwise:
            I_net_id(find(I_net_id==intersect(I_net_id,I_id_net)))=[];
            % get station info 
            sta_info.name = id_ref; 
            sta_info.lon  = sta.lon(ksta);
            sta_info.lat  = sta.lat(ksta); 
            sta_info.elev = sta.elev(ksta);
        end
        function sta=get_station_list(h)
            [sta_ Ista]=unique(h.id(:));
            sta.lon=double(h.lon(Ista));
            sta.lat=double(h.lat(Ista));
            sta.elev=double(h.elev(Ista));
            sta.sta=h.id(Ista);
        end

        %----------------------------------------------
        %        read correlations from their indice  
        %----------------------------------------------
        function [cc_matrix ref]=read_cc_matrix(h,kpath,p1,p2,kcmp)
        % read all daily correlation for a single station pair 
            if ~exist('kcmp') ; kcmp=1; end
            [filename I]=h.get_filename_from_indice(kpath);
            path_ = h.get_path_within_h5(kpath,kcmp);
            nwf =numel(h.t)     ;
            nday=numel(h.date1) ;
            cc_matrix=h5read(filename,['/cc',path_]); %1,I,1,kcmp],[nwf,1,nday,1]);
            cc_matrix=squeeze(double(cc_matrix));
            if nday>5000 
                cc_matrix=s2d.filter_low_memory(cc_matrix,p1,p2,h.tau);
            else
                cc_matrix=s2d.filter(cc_matrix,p1,p2,h.tau);
            end
            ref=h.read_cc(kpath,p1,p2,'cmp',kcmp) ;
        end
        
        function cc=read_cc_sta_vs_others(h,ksta,p1,p2,varargin)
            in.cmp=1     ; % vector of components to be read
            in.day=false ; % vector of date to be read (false = reference)
            in.taper=200 ; % taper length (npts) applied at both side of each cc 
            in=lang.parse_options(in,varargin);
            
            sta=h.get_station_list; 
            id_ref=sta.sta{ksta}  ;
            %get all correlation id_ref -> network : 
            I_id_net=find(strcmp(h.id(1,:),id_ref)==1); 
            cc_id_net=h.read_cc(I_id_net,p1,p2,in)    ;
            %get all correlation network -> id_ref : 
            I_net_id    =find(strcmp(h.id(2,:),id_ref)==1); 
            %   we need to remove the autocorr since we would get it twice otherwise:
            I_net_id(find(I_net_id==intersect(I_net_id,I_id_net)))=[];
            cc_net_id   =h.read_cc(I_net_id,p1,p2,in);
            cc_net_id.cc= flipud(cc_net_id.cc);
            %merge the two matrices : 
            cc=struct();
            cc.cc   =[cc_id_net.cc   cc_net_id.cc]      ;
            cc.dist =[cc_id_net.dist cc_net_id.dist]    ;
            cc.titre=[cc_id_net.titre ; cc_net_id.titre];
            cc.p1   = cc_id_net.p1 ; 
            cc.p2   = cc_id_net.p2 ;
            %sort the result by distance : 
            [~, J]   = sort(cc.dist); 
            cc.cc   = cc.cc(:,J)    ; 
            cc.dist = cc.dist(J)    ; 
            cc.titre= cc.titre(J)   ; 
        end	
        
        function cc=read_cc(h,kpath,p1,p2,varargin)
            in.cmp=1     ; % vector of components to be read
            in.day=false ; % vector of date to be read (false = reference)
            in.taper=200 ; % taper length (npts) applied at both side of each cc 
            in=lang.parse_options(in,varargin);
            
            npath=numel(kpath) ;
            nper =numel(p1);
            ncmp =numel(in.cmp);
            nday =numel(in.day);
            nwf  =numel(h.t);
            ntr  =npath*nper*ncmp*nday;

            cc      =struct();
            cc.cc   =zeros(nwf,ntr);
            cc.titre=cell(ntr,1) ;
            cc.dist =zeros(1,ntr);
            cc.p1   =zeros(1,nper);
            cc.p2   =zeros(1,nper);
            cc.ref_nstack=zeros(1,ntr);
            cc_az   = zeros(1,ntr); 
            kcc=1;	
            
            for iper=1:nper
                cp1=p1(iper);
                cp2=p2(iper);
                for ipath = kpath 
                    for icmp = in.cmp 
                        for iday = in.day 
                            cc.cc(:,kcc)   =h.read_single_cc(ipath,cp1,cp2,icmp,iday,in.taper);
                            cc.titre{kcc}=h.get_title(ipath,cp1,cp2,icmp,iday);
                            cc.dist(kcc) =h.dist(ipath);
                            cc.p1(kcc)   =cp1; 
                            cc.p2(kcc)   =cp2; 
                            cc.ref_nstack(kcc)=h.ref_nstack(ipath);
                            cc.az(kcc)   = h.az(ipath);
                            kcc=kcc+1;
                        end
                    end
                end
            end
        end	
        
        
        function cc=read_single_cc(h,kpath,p1,p2,kcmp,kday,ntaper)
        % read a single correlation (daily or ref is day is not set)
            ref=false;
            if ~exist('kday') | isempty(kday);	
                ref=true ; 
            elseif kday == false  ;
                ref=true ; 
            end
            if ~exist('kcmp') | isempty(kcmp) 
                kcmp=1   ;
            end
            if ~exist('ntaper') ; 
                ntaper=0; 
            end
            
            [filename I]=h.get_filename_from_indice(kpath);
            path_ = h.get_path_within_h5(kpath,kcmp);
            if ref 
                cc=h5read(filename,['/ref',path_]); %,[1,I,kcmp],[numel(h.t),1,1]);
            else 
                cc=h5read(filename,['/cc',path_],[1,kday],[numel(h.t),1]);
            end
            cc=double(cc);
            if exist('p1')
                if ntaper >0 
                    cc=s2d.taper_edges(cc,ntaper);
                end		 
                cc=s2d.filter(cc,p1,p2,h.tau,4);
            end
        end
        
        function path_=get_path_within_h5(h,I,kcmp)
            path_=['/',h.id{1,I},'/',h.id{2,I},'/',h.cmp{kcmp}];
        end
        
        function titre=get_title(h,kpath,p1,p2,kcmp,kday) 
            titre  =[h.id{1,kpath},'-',h.id{2,kpath}];
            titre  =[titre,':   ',str.num2str(round(h.dist(kpath)),3),' km'];
            if exist('p1');
                titre  =[titre,'   ',str.num2str(p1,2),'-',str.num2str(p2,2),'s'];
            end
            if exist('kcmp')
                titre  =[titre,'   ',h.in.cc_cmp{kcmp}];
            end
            if exist('kday')
                if kday==false 
                    titre =[titre,'   ref'];
                else
                    titre =[titre,'  ' datestr(h.date1(kday),'yyyy/mm/dd HH:MM')];
                    titre =[titre,'  ' datestr(h.date2(kday),'yyyy/mm/dd HH:MM')];
                end
            end
        end	  	
        
        function [filename I_within_file]=get_filename_from_indice(h,kpath)
        % return the filename of the correlation #I 
            I_file=find(h.file_I(2,:) >= kpath,1,'first') ;
            filename=h.file{I_file};
            I_within_file = kpath-h.file_I(1,I_file)+1;
        end 
        
        %----------------------------------------------
        %                 constructor                  
        %----------------------------------------------
        function reformat_filename(h,in) 
        % first reconstruct the name of the dir where the xcorr files are
            if in(end)=='5'                 % we read a pydb*h5 files from pycorr.live.disp_ for instance
                aa = strsplit(in,'/') ;     % remove the name of the .h5 to keep only the directory name
                in = join(aa(1:end-1),'/');
                in = in{1};
            end
            if in(end)=='/';                 % the user  specified a directory name with a trailing /
                in = in(1:end-1);           % so we removed it.               
            end
            nfile = numel(h.file);    
            for ifile =1 :nfile
                aa=strsplit(h.file{ifile},'/') ;
                h.file{ifile}=[in,'/',aa{end}];
            end
        end
        %--------------------
        function h=live_md(in)				
            is_dir=isdir(in);
            if is_dir== 1  %working with a C1 dir :
                if ~exist([in,'/db.mat']);
                    dispc('we need to build the db.mat file','c','b');
                    db=cmd_build_db(in);
                    try  % in case we work in s.o else directory and have no write access 
                        save([in,'/db.mat'],'-struct','db');
                    catch	
                        save(['db.mat'],'-struct','db');
                    end
                end
                if exist('db')~=1
                    h=struct_.copy_field_to(load([in,'/db.mat']),h);
                else
                    h=struct_.copy_field_to(db,h);
                end
            else isstr(in) ;
                if strcmp(in(end-2:end),'.h5') % this is pydb_.h5 file : do nothing ?
                    dispc('  pycorr.live_md : read correlations metadata','c','n')
                    md=read_md_from_pydb_h5_file(in);
                    
                    if isfield(md,'date1_ordinal')
                        md.date1 = md.date1_ordinal; 
                        md.date2 = md.date2_ordinal; 
                        md = rmfield(md,'date1_ordinal'); 
                        md = rmfield(md,'date2_ordinal'); 
                    end
                    h=struct_.copy_field_to(md,h);
                elseif strcmp(in(end-2:end),'mat') % this is db_disp*mat file :
                                                   %md = read_md_from_db_disp
                                                   %h=struct_.add_fields_except(matfile(in),{'disp','Properties'},h);
                end
                h.reformat_filename(in) % in case we run pycorr.live from another directory than the 1st time when we
                                      % constructed the db.mat file.
            end
        end
    end
end


%--------------------------------------------------------------------------------------
%
% the following functions are outside the class, since they do manipulate h5 files and 
% not the class properties (they do not access h.*)
%----------------------------------------------------------------------------------------

function md=read_md_from_pydb_h5_file(filename)
	%this function is called when reading a pydb_dvv_*.h5 file. In this case 
	%the h5 file already contains dist/az/baz/file/file_I fields

	%read pydb_h5 file in the same way we extract md from xcorr*h5 file :
	md=cmd_read_metadata(filename)	    ;
	%do the same cleaning as when we concatenate xcorr metadata :
	if isfield(md,'date_1_ordinal')
		md.date1 = md.date1_ordinal+366    ;
		md.date2 = md.date2_ordinal+366     ; 
		md = rmfield(md,'date1_ordinal')    ;
		md = rmfield(md,'date2_ordinal')    ;
	end
	%add the file and file_I information : 
	md.file   = h5read(filename,'/file');
	md.file_I = h5read(filename,'/file_I')+1;
	md.id=cellfun(@celltrim,md.id)       ;
    end

function md=cmd_build_db(in)
    h5_list=os.dir([in,'/[xv][cs]*h5']);
    %put each h5 file metadata into md_all array : 
    for ih5 = h5_list
        dispc(ih5,'c','b');
        if ~exist('md_all') 
            md_all=cmd_read_metadata(ih5{1})	;
        else
            md_all(end+1)=cmd_read_metadata(ih5{1});
        end
    end
    %flatten md_all : 
    md=md_all(1)         ;
    md.id  =[md_all.id]  ;
    md.lat =[md_all.lat] ;
    md.lon =[md_all.lon] ;
    md.elev=[md_all.elev];		
    if isfield(md,'ref_nstack') 
        md.ref_nstack=[md_all.ref_nstack];
    end
    %compute the dist, az and baz of each station pair : 
    npath = size(md.lon,2);
    md.dist=zeros(1,npath);
    md.az  =zeros(1,npath);
    md.baz =zeros(1,npath);
    for ipath =1 : npath 
        lon1=md.lon(1,ipath); 
        lon2=md.lon(2,ipath);
        lat1=md.lat(1,ipath);
        lat2=md.lat(2,ipath);
        [md.dist(ipath) md.az(ipath) md.baz(ipath)]=m_idist(lon1,lat1,lon2,lat2);
    end
    md.dist=md.dist/1000; 
    md.dist(find(isnan(md.dist)==1))=0;
    %store file name : 
    nfile    =numel(md_all) ;
    md.file  =cell(1,nfile) ;
    md.file_I=zeros(2,nfile);
    I1=1;
    for ifile =1 :nfile 
        md.file(ifile)=h5_list(ifile)    ;
        I2=I1+size(md_all(ifile).lon,2)-1;
        md.file_I(:,ifile)=[I1 I2]       ;
        I1=I2+1                          ;
    end
    %convert python date to matlab dates : 
    if isfield(md,'date1_ordinal')        ; %does not exist for c3 conveted into c1
        md.date1 = md.date1_ordinal+366     ;
        md.date2 = md.date2_ordinal+366     ; 
        md = rmfield(md,'date1_ordinal')    ;
        md = rmfield(md,'date2_ordinal')    ;
    end
    %remove white space at the end of ID name : 
    md.id=cellfun(@celltrim,md.id)      ;
end


function md=cmd_read_metadata(h5_file)
    md   = h5.read_group(h5_file,'/md') ;
    md.in= h5.read_group(h5_file,'/in_');
    md   = struct_.copy_field_to(h5.read_group(h5_file,'/md_c'),md);
    try 
        md.ref_nstack = h5read(h5_file,'/ref_nstack')';
    catch
        dispc(['  ',h5_file,' : ref_nstack was not found'],'c','n')
    end
    if isfield(md,'I1')     ; %does not exist for c3 converted into c1
        md.in.I1=md.I1      ;
        md.in.I2=md.I2      ;
        md =rmfield(md,'I1');
        md =rmfield(md,'I2');
    end
end	


function str_=celltrim(str_)
	str_={str_(find(str_~=char(0)))};
end
