%----------------------------------------------------------------------------
% class to handle Cvideo files that contains vector of dvv per station.
% These files are typically outputted by : 
%    h=pycorr.dvv_(cpydb);
%    h.dvv_build_video_around_station(in) 
%    => call dvv_.dvv_p_around_station(ista) in a loop on stations
%    it take all dvv measurements performed around each station pre-process them using dvv_selection and
%    averaged them. The resulting dv/v per station are saved in Cvideo/*mat files. 
%-----------------------------------

classdef dvv_video < handle 
    properties
        dvv      =[] ; % [218×278] les dv/v
        coh      =[] ; % [218×278] coherence de chaque mesire
        dvv_std  =[] ; % [218×278] std des mesures des N trajets moyennees autour de chaque station 
        rms      =[] ; % [218x278] incertitude sur les dv/v d'apres weaver et al.m
        date2    =[] ; % [218×1 double]
        date1    =[] ;
        dist_area=[] ; % [1×278 double]
        npath =[]    ; % [1×278 double]
        fname =[]    ;  
        in=struct()  ;
        sta=struct() ;
    end
    methods
        function p_dist_area_gmap(h,varargin) 
            in.ss=80; 
            in.width=0.7750;
            in.caxis=[0 300];
            in.alpha=0.4;
            in.type='terrain';
            in.resize=2; 
            in.scale=2;
            in.refresh=true;
            in=lang.parse_options(in,varargin);
            scatter(h.sta.lon,h.sta.lat,in.ss,h.dist_area,'fill');
            
            %modify the width of the figure useful for gmap: 
            aa=get(gca,'Position');
            dx=aa(3) - in.width;
            set(gca,'Position',[aa(1)+dx/2 aa(2) in.width aa(4)])
            
            gmap.map([],[],'type',in.type,'alpha',in.alpha,'scale',in.scale,'resize',in.resize); 
            caxis(in.caxis)         ;
            hc=colorbar             ;  
            colorTitleHandle = get(hc,'Title');
            set(colorTitleHandle ,'String','km');
            
            p.cb_make_me_thinner(hc);
            titre='spatial averaging done around each station [km]';
            title(titre,'fontweight','bold','fontsize',16);
            colormap(flipud(jet));
        end
        function ha=p_map_and_dvv(h,cdate,ksta,varargin);
            in.ylim=[-1 1]*1e-3;
            in.rms=true;
            in.col='k';
            in.lw = 2 ;
            in.rms_col = [0.8 0.8 1];
            in.rms_edge_col=[0 0 0] ;
            in.rms_alpha=[0.2]      ;
            in.eq_date=[]           ;
            in.caxis = [-1 1]*1e-3  ;
            in.eq_lon =[]           ;
            in.eq_lat =[]           ;
            in.eq_catalog =[]       ;
            in.gmap=false           ; 
            in.gm_width=0.5         ; % if we use gmap : width of the plot
            in.gm_alpha=0.5           ; % if we use gmap : transparency of the plot
            in.gm_type='terrain'    ;
            in.refresh= true        ;
            in.print =false         ;
            in.print_video=[]       ; %empty for false
            in.lon=[]               ;
            in.lat=[]               ;
            in.dpi=150              ;
            in.func=[]              ;
            in = lang.parse_options(in,varargin);

            if in.print
                aa=strsplit(h.fname,'/');
                db=dbstack;
                if numel(db)==2 % called by another script
                    out_dir =['plot/',db(2).name,'__dvv_map_at_date/',aa{2}(1:end-4)];
                else
                    out_dir = ['plot/dvv_map_at_date/',aa{2}(1:end-4)];
                end
                os.mkdir(out_dir)
                if numel(in.print_video)>0;
                    out_file=[out_dir,'/',str.num2str(in.print_video,5),'__video.png'];
                else
                    out_file = [out_dir,'/',datestr(cdate,'yyyy_mm_dd'),'.png'];
                end
                %p.png('big','dpi',in.dpi,'out_file',out_file);
                if exist(out_file,'file');
                    ha=false;
                    dispc([out_file,' already exists'],'r','r')
                    return
                end
            end

            if in.refresh
                ha= p.split_map_trace;
            else
                ha=in.ha;
            end
            subplot(ha(1)); 
            if ~in.gmap
                I=h.p_map_at_date(cdate,'black',true,'caxis',in.caxis);
            else
                I=h.p_gmap_at_date(cdate,'caxis',in.caxis,'width',in.gm_width,'alpha',in.gm_alpha,'type',in.gm_type,'refresh',in.refresh,'lon',in.lon,'lat',in.lat);
            end
            if numel(in.eq_catalog) > 0
                I_eq = find(in.eq_catalog.date >= h.date1(I) & in.eq_catalog.date <= h.date2(I));
                plot(in.eq_catalog.lon(I_eq),in.eq_catalog.lat(I_eq),'p','color','k', 'MarkerFaceColor','y')
                I_eq = find(in.eq_catalog.date >= h.date1(I) & in.eq_catalog.date <= h.date2(I) & in.eq_catalog.Mw >5);
                plot(in.eq_catalog.lon(I_eq),in.eq_catalog.lat(I_eq),'p','color','k', 'MarkerFaceColor','r','markersize',12)
            
            end
            %map.station(h.sta.lon(ksta),h.sta.lat(ksta),'mark','^','msize',10)
            if numel(in.eq_lon)>0 
                keyboard
            end
            
            subplot(ha(2))
            cla;
            h.p_dvv_at_station(ksta,in);
            ylim(in.caxis);
            ylimit=ylim;
            a=[h.date1(I) h.date2(I) h.date2(I) h.date1(I)];
            b=[ylimit(1) ylimit(1) ylimit(2) ylimit(2)];
            patch(a,b,'r','facealpha',0.2);
            
            %optionnaly exectute an external function :
            if ~isempty(in.func)
               in.func(ha)
            end
            
            if in.print
                %aa=strsplit(h.fname,'/');
                %db=dbstack;
                %if numel(db)==2 % called by another script
                %    out_dir =['plot/',db(2).name,'__dvv_map_at_date/',aa{2}(1:end-4)];
                %else
                %    out_dir = ['plot/dvv_map_at_date/',aa{2}(1:end-4)];
                %end
                %os.mkdir(out_dir)
                %if numel(in.print_video)>0;
                %    out_file=[out_dir,'/',str.num2str(in.print_video,5),'__video.png'];
                %else
                %    out_file = [out_dir,'/',datestr(cdate,'yyyy_mm_dd'),'.png'];
                %end
                p.png('big','dpi',in.dpi,'out_file',out_file);
            end
        end
        function p_dvv_at_station(h,ksta,varargin);
            in.title=true;
            in.ylim=false;
            in.rms=false ; 
            in.col='k';
            in.rms_col='k';
            in.rms_edge_col='w';
            in.lw=1;
            in.rms_alpha=0.2;
            in.eq_date = [];
            in= lang.parse_options(in,varargin);

            mk_title=in.title; 
            in = rmfield(in,'title');
            if in.rms 
                in.rms = h.rms(:,ksta);
            end
            p.date_plot(h.date2,h.dvv(:,ksta),in);%
            set(gca,'fontsize',12); 
            if in.eq_date
                p.date_add_line_at_date(in.eq_date);
            end
            if mk_title
                titre =['dv/v at ',h.sta.sta{ksta}];
                title(titre,'fontweight','bold','fontsize',16);
            end
            
        end
        function [I hc]=p_gmap_at_date(h,date_,varargin);
            in.ss=80; 
            in.caxis=[-1 1]*1e-3;
            in.cmap='jet';
            in.width=0.7750; 
            in.alpha=0.4;
            in.type='terrain';
            in.resize=2; 
            in.scale=2;
            in.refresh=true;
            in.lon=[];
            in.lat=[]; 
            in = lang.parse_options(in,varargin);

            if numel(in.lon)>1 
                I_sta=h.get_indice_of_station_within_coordinate(in.lon,in.lat);
            else
                I_sta=[1:numel(h.sta.lon)];
            end
            I=find(h.date2>=date_,1,'first')                      ;
            scatter(h.sta.lon(I_sta),h.sta.lat(I_sta),in.ss,h.dvv(I,I_sta),'fill');
            aa=get(gca,'Position');
            dx=aa(3) - in.width;
            set(gca,'Position',[aa(1)+dx/2 aa(2) in.width aa(4)])
            if in.refresh
                gmap.map([],[],'type',in.type,'alpha',in.alpha,'scale',in.scale,'resize',in.resize); 
            end
            caxis(in.caxis)         ;
            hc=colorbar             ;  
            p.cb_make_me_thinner(hc);
            %titre=['dvv : ',datestr(h.date2(I),'yyyy-mm-dd')];
            titre=[datestr(h.date1(I),'yyyy-mm-dd'),' - ',datestr(h.date2(I),'yyyy-mm-dd')];
            title(titre,'fontweight','bold','fontsize',16);
            colormap(in.cmap);

            
        end
        function I=p_map_at_date(h,date_,varargin)
            in.ss=80; 
            in.caxis=[-1 1]*1e-3;
            in.black=false;
            in.cmap='jet';
            in = lang.parse_options(in,varargin);
            
            h.p_map('topo',false,'p_sta',false,'black',in.black)  ;
            I=find(h.date2>=date_,1,'first')                      ;
            m_scatter(h.sta.lon,h.sta.lat,in.ss,h.dvv(I,:),'fill');
            caxis(in.caxis)         ;
            hc=colorbar             ;  
            p.cb_make_me_thinner(hc);
            titre=['dvv : ',datestr(h.date2(I),'yyyy-mm-dd')];
            %titre=['dvv : ',datestr(h.date1(I),'yyyy-mm-dd'),'-',datestr(h.date2(I),'yyyy-mm-dd')];
            title(titre,'fontweight','bold','fontsize',16);
            colormap(in.cmap);
        end
        function p_map(h,varargin)
            in.dx=0.1;
            in.dy=0.1; 
            in.topo=false;
            in.p_sta= true; 
            in.sta_size=16;
            in.sta_sym='^';
            in.black=false ;
            in = lang.parse_options(in,varargin);
            
            inm.topo=in.topo;
            if in.black 
                inm.coastcol=[1 1 1]; 
                inm.mgrid_color=[1 1 1]*0.5;
            else
                inm.coastcol=[0 0 0]; 
            end
            
            hold on;
            lon=[min(h.sta.lon)-in.dx max(h.sta.lon)+in.dx];
            lat=[min(h.sta.lat)-in.dx max(h.sta.lat)+in.dx];
            map.background(lon,lat,inm);
            if in.p_sta
                map.station(h.sta.lon,h.sta.lat,'mark',in.sta_sym,'color','w','msize',in.sta_size);
            end
        end
        function dvv_correct_offset(h,date1,date2,varargin)
            I=find(h.date2>=date1 & h.date2<=date2);
            offset = mean(h.dvv(I,:),'omitnan'); 
            Inan=find(isnan(offset)==1);
            dispc(['  ',num2str(numel(Inan)),' paths set to nan when correcting offset'],'c','n');
            h.dvv = h.dvv-offset;
        end
        function I_sta=get_indice_of_station_within_coordinate(h,lon,lat)
            I1=find(h.sta.lon >= lon(1) & h.sta.lon <= lon(2));
            I2=find(h.sta.lat >= lat(1) & h.sta.lat <= lat(2));
            I_sta = intersect(I1,I2);
        end
        function h=dvv_video(fname);
            out=load(fname)   ;   
            h.dvv = out.dvv   ;
            h.coh = out.coh   ;
            h.rms = out.rms   ;
            h.dvv_std = out.dvv_std ;
            h.date2   = out.date2   ;
            h.date1   = out.date1   ; 
            h.dist_area = out.dist_area; 
            h.npath     = out.npath    ;
            h.in        = out.in       ;
            h.fname     = fname        ;
            h.sta       = out.sta      ;
        end
    end
end