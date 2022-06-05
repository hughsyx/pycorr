% class to manipulate a subset of dvv measurement. Here we loose completely the link to the original correlation. 
% if we add more components => put them on the third dimension
classdef dvv_selection < handle
    properties
        dvv  =[]   ;% [202×293 double]
        coh  =[]   ;% [202×293 double]
        lon  =[]   ;% [2×293 double]
        lat  =[]   ;% [2×293 double]
        id   =[]   ;% {2×293 cell}
        dist =[]   ;% [1×293 single]
        az   =[]   ;% [1×293 single]
        baz  =[]   ;% [1×293 single]
        date2 =  [];% 
        date1=   [];%
        fname=''   ;
        in=struct()     ; % input parameters used to compute dv/v 
        sta=struct()    ; % list of stations kept 
        sta_ori=struct(); % original list of stations when calling dvv_selection
    end
    methods % method used to plot the measurements
        function [I_date]=plot_gmap_of_path_at_date(h,cdate,varargin); 
            in.dx=0.1;
            in.dy=0.2;
            in.caxis=[-1 1]*1e-3;
            in.lw =1;
            in.width=[];
            in = lang.parse_options(in,varargin);
            
            % prepare the plot :
            I_date= find(h.date2>=cdate,1,'first');
            hc=p.line_with_color(h.lon,h.lat,h.dvv(I_date,:),in.caxis,'lw',in.lw);
            p.cb_make_me_thinner(hc);
            if ~isempty(in.width);
                aa=get(gca,'Position');
                width = aa(3)*in.width;
                dx = aa(3)-width;
                set(gca,'Position',[aa(1)+dx/2 aa(2) width aa(4)]);
            end
            gmap.map([],[],'type','hybrid','alpha',0.5,'scale',2,'resize',2) 
            %h.p_gmap('msize',8)
            %title : 
            titre = [datestr(cdate,'yyyy/mm/dd'),' dv/v per path'];
            title(titre,'fontsize',14,'fontweight','bold')
            
        end
        function [I_date]=plot_map_of_path_at_date(h,cdate,varargin);
            in.black=true; 
            in.dx=0.1;
            in.dy=0.2;
            in.caxis=[-1 1]*1e-3;
            in.lw =1;
            in = lang.parse_options(in,varargin);
            
            %plot map background :
            lon=[min(h.sta.lon)-in.dx max(h.sta.lon)+in.dx];
            lat=[min(h.sta.lat)-in.dx max(h.sta.lat)+in.dx];
            if in.black
                map.background(lon,lat,'topo',false,'coastcol',[0.5 0.5 0.5],'mgrid_color',[0.5 0.5 0.5]);
            else
                map.background(lon,lat,'topo',false);
            end
            % prepare the plot :
            I_date= find(h.date2>=cdate,1,'first');
            map.line_with_color(h.lon,h.lat,h.dvv(I_date,:),in.caxis,'lw',in.lw);
            map.station(h.sta.lon,h.sta.lat);            
            %title : 
            titre = [datestr(cdate,'yyyy/mm/dd'),' dv/v per path'];
            title(titre,'fontsize',14,'fontweight','bold')
        end
        function [ha] = plot_path_map_and_dvv(h,cdate,varargin)
            in.black=true; 
            in.dx=0.1;
            in.dy=0.2;
            in.caxis=[-1 1]*1e-3;
            in.lw =1;
            in.ylim=[-1 1]*1e-3;
            in.eq_date=false;
            in = lang.parse_options(in,varargin);
            
            ha=p.split_map_trace;
            subplot(ha(1))
            [I]=h.plot_map_of_path_at_date(cdate,in);
            
            subplot(ha(2))
            hold on
            p.date_plot(h.date2,mean(h.dvv,2,'omitnan'),'ylim',in.ylim);
            
            %rms=h.compute_rms;
            if in.eq_date
                p.date_add_line_at_date(in.eq_date);
            end
            % add a box that tells the date:
            ylimit=ylim;
            a=[h.date1(I) h.date2(I) h.date2(I) h.date1(I)];
            b=[ylimit(1) ylimit(1) ylimit(2) ylimit(2)];
            patch(a,b,'c','facealpha',0.2);
            
        end
        function [ha]=plot_map_and_matrix(h,varargin)
            in.xlim=false                     ; 
            in.datetick_fmt='yyyy/mm/dd'      ;
            in.p_mean_ylim=[-6 2]*1e-4        ;
            in.print=false                    ;
            in.print_name=''                  ;
            in.title=''                       ;
            in.eq_date=[]                     ;
            in.black=false                    ;
            in.gmap=false                     ; 
            in.gwidth=1                       ;
            in.lw=1;                      
            in=lang.parse_options(in,varargin);
            
            clf; 
            ha=p.split_map_nsubplot(5)        ;
            %selected station map :
            subplot(ha(1));
            if in.gmap 
                h.p_gmap('width',in.gwidth);
            else
                h.p_map;
            end
            if numel(in.title)>0
                titre=in.title;
            else
                titre=h.dvv_get_title_generic;
            end
            title(titre,'fontsize',14,'fontweight','bold');
            % mean dvv :
            subplot(ha(2))
            h.p_mean_dvv('xlim',in.xlim,'datetick_fmt',in.datetick_fmt,'ylim',in.p_mean_ylim,'p_med',false,'lw',in.lw);
            
            %subplot(ha(3))
            %std_ = std(h.dvv,[],2,'omitnan');
            %p.date_plot(h.date2,std_,'ylim',[0 1]*1e-3,'xlim',in.xlim,'date_fmt',in.datetick_fmt);
            %title('std deviation of dv/v measurements','fontweight','bold')
            subplot(ha(3))
            p.date_matrix(h.date2,h.dist,h.dvv','caxis',[-1 1]*1e-3);
            ylabel('path sorted by distance with imagesc','fontweight','bold')
            title('dvv : matrix vs distance','fontweight','bold'); 
            
            % proba plot :
            subplot(ha(4))
            %p.date_matrix_proba(h.date2,h.dvv,h.in.stretch(1:2:end));
            title('probability of measurement','fontweight','bold');
            % mean coherency
            subplot(ha(5))
            h.p_mean_coh('xlim',in.xlim,'datetick_fmt',in.datetick_fmt,'ylim',[0 1],'lw',in.lw);
            
            subplot(ha(6));
            p.date_matrix(h.date2,h.dist,h.coh','caxis',[0 1]);
            ylabel('path sorted by distance with imagesc','fontweight','bold')
            title('coherency matrix vs distance','fontweight','bold');
            if numel(in.eq_date)>0 
                for iha=2:6
                    subplot(ha(iha))
                    p.date_add_line_at_date(in.eq_date);
                end
            end
            if in.black==true  & in.gmap==false
                subplot(ha(1));
                m_grid('color','white');
            end
            
            if in.print
                aa=strsplit(h.fname,'/');
                outdir = ['plot/p_around_station/',aa{2}(1:end-3)];
                os.mkdir(outdir);
                out_file =[outdir,'/',in.print_name,'.png'];
                p.png('big','dpi',100,'out_file',out_file);
            end
        end
        function p_mean_coh(h,varargin)
            in.xlim=false; 
            in.ylim=[0 1] ;
            in.datetick_fmt='yyyy/mm/dd';
            in.title=true;
            in.lw=1;
            in = lang.parse_options(in,varargin);
            
            plot(h.date2,mean(h.coh','omitnan'),'linewidth',in.lw);
            
            ylim(in.ylim) ;

            if in.xlim==false
                xlim([h.date2(1) h.date2(end)]);
            else
                xlim(in.xlim);
            end
            datetick('x',in.datetick_fmt,'keeplimits'); 
            grid minor;
            if in.title
                title('mean coda coherency','fontweight','bold');
            end
        end
        function p_mean_dvv(h,varargin);
            in.xlim=false; 
            in.ylim=double([h.in.stretch(1) h.in.stretch(end)]);
            in.datetick_fmt='yyyy/mm/dd';
            in.title=true;
            in.p_med=false;
            in.lw=1;
            in = lang.parse_options(in,varargin);

            plot(h.date2,mean(h.dvv','omitnan'),'linewidth',in.lw);
            if in.p_med
                hold on;
                plot(h.date2,median(h.dvv','omitnan'),'--','linewidth',in.lw);
                %legend('mean','median')
            end
            if numel(in.ylim)==2
                ylim(in.ylim) ;
            else
                try
                    mdvv=mean(h.dvv','omitnan');
                    ylim([min(mdvv) max(mdvv)]);
                catch
                    ylim([-5 5]*1e-4);
                end
            end
            if in.xlim==false
                xlim([h.date2(1) h.date2(end)]);
            else
                xlim(in.xlim);
            end
            datetick('x',in.datetick_fmt,'keeplimits'); 
            grid minor;
            if in.title
                title('mean dvv','fontweight','bold');
            end
        end
        function p_gmap(h,varargin)
            in.xlim=[];
            in.ylim=[];
            in.width=[];
            in.msize=12;
            in.p_sta_ori=true;
            in = lang.parse_options(in,varargin);
            hold on
            if in.p_sta_ori 
                plot(h.sta_ori.lon,h.sta_ori.lat,'k^','MarkerFaceColor',[0.5 0.5 0.5],'Markersize',in.msize) 
            end
            plot(h.sta.lon,h.sta.lat,'k^','MarkerFaceColor','r','markersize',in.msize) 
            
            if ~isempty(in.width);
                aa=get(gca,'Position');
                width = aa(3)*in.width;
                dx = aa(3)-width;
                set(gca,'Position',[aa(1)+dx/2 aa(2) width aa(4)]);
            end
            if ~isempty(in.xlim);
                xlim(in.xlim);
                ylim(in.ylim);
            end
            gmap.map([],[],'type','terrain','alpha',1,'scale',2,'resize',2) 
        end
        function p_map(h,varargin)
            in.dx=0.1;
            in.dy=0.1;
            in.p_sta_ori=false;
            in.topo=true;
            in.black=false;
            in = lang.parse_options(in,varargin);
            hold on;
            
            if in.p_sta_ori
                lon=[min(h.sta_ori.lon)-in.dx max(h.sta_ori.lon)+in.dx];
                lat=[min(h.sta_ori.lat)-in.dx max(h.sta_ori.lat)+in.dx];
            else
                lon=[min(h.sta.lon)-in.dx max(h.sta.lon)+in.dx];
                lat=[min(h.sta.lat)-in.dx max(h.sta.lat)+in.dx];
            end
            map.background(lon,lat,'topo',in.topo,'black',in.black);
            if in.p_sta_ori
                map.station(h.sta_ori.lon,h.sta_ori.lat,'mark','^','color','w');
            end
            map.station(h.sta.lon,h.sta.lat,'mark','^','color','r');
        end
        function titre=dvv_get_title_generic(h)
            titre=[];
            titre = ['dv/v ',num2str(h.in.p1),'-',num2str(h.in.p2),'s'];
            titre = [titre,'   stack ',num2str(h.in.stack),'-',num2str(h.in.stack_shift)];
            titre = [titre,'   coda length ',[num2str(h.in.tf),'s']];
        end

    end
    methods %used to construct the object and to select path and measurements:
            %function derivate_dvv(h,varargin)
            %npath = numel(h.dist);
            %ipath=1
            %keyboard
            %dvv = nan(size(h.dvv));
            %for ipath=1:npath
            %    toto=fillmissing(h.dvv(:,ipath),'linear');
            %    dvv(:,ipath)=s1d.derivate(toto,h.date2);
            %end
            %end
        function rms=compute_rms(h)
            wc=1./mean(double([h.in.p1 h.in.p2])); %Hz
            T= 1./(double(h.in.p2)-double(h.in.p1))      ; %T [Hz]
            dist_mean=mean(h.dist);
            t1_mean=dist_mean/h.in.v1+h.in.coda_dp*h.in.p2;
            t2_mean=t1_mean+h.in.tf;%h.mt.it(1).tf;
            ncmp = size(h.coh,3);
            ndate = size(h.coh,1);
            rms = nan(ndate,ncmp);
            for icmp =1 : ncmp 
                rms1= sqrt((6*T*sqrt(pi/2)))/(wc^2*(t2_mean^3-t1_mean^3))  ;
                cc_mean= mean(squeeze(h.coh(:,:,icmp))','omitnan');
                rms2= (sqrt(1-cc_mean.^2))./(2*cc_mean);
                nmeas=sum(~isnan(h.dvv(:,:,icmp)),2);
                %rms2= (sqrt(1-x.cc_mean(:,icmp).^2))./(2*x.cc_mean(:,icmp));
                rms(:,icmp)=rms1*rms2'./sqrt(nmeas);
            end
        end
        function keep_measurement_with_cc_btw(h,cc)
            I=find(h.coh(:)< cc(1) | h.coh(:) > cc(2));
            h.dvv(I)=nan;
            h.coh(I)=nan;
            cc_str=[num2str(cc(1)),'-',num2str(cc(2))];
            msg1=['  setting all dv/v measurements with a cc outside ',cc_str,' to nan'];
            msg2=['  => ',num2str(numel(I)),' meas out of ',num2str(numel(h.dvv(:))),' set to nan'];
            dispc(msg1,'c','b');
            dispc(msg2,'c','n');
        end
        function rm_path_with_av_coh_inside(h,cc,varargin)
            in=struct(); 
            in.date1=datenum(2000,1,1); 
            in.date2=datenum(2050,1,1); 
            in=lang.parse_options(in,varargin);
            I_date = find(h.date2 >= in.date1 & h.date2 <= in.date2); 
            coh_mean=mean(h.coh(I_date,:),'omitnan');
            
            I_path = find(coh_mean >= cc(1) & coh_mean <= cc(2));
            h.rm_path(I_path);
        end
        function  rm_path_with_less_than_nday_of_data(h,nday)
            % remove path that have less than nday of measurements (we correct the numbef of meas
            % by stackshift so that it is in day
            nday_per_path = numel(h.date2)-sum(isnan(h.dvv));
            nday_per_path = nday_per_path*double(h.in.stack_shift);
            I_path=find(nday_per_path<nday);
            h.rm_path(I_path);
            msg1=['  removing all path with less than ',num2str(nday),' days of data'];
            msg2=['  => ',num2str(numel(I_path)),' paths removed out of ',num2str(numel(nday_per_path))];
            %dispc('  --','c','n');
            dispc(msg1,'c','b');
            dispc(msg2,'c','n');
        end
        function rm_path(h,I_path)
            h.dvv(:,I_path)=[];
            h.coh(:,I_path)=[];
            h.lon(:,I_path)=[];
            h.lat(:,I_path)=[];
            h.id(:,I_path)=[];
            h.dist(I_path)=[];
            h.az(I_path)=[];
            h.baz(I_path)=[];
        end
        function correct_each_dvv_from_mean(h,date1,date2)
            fmt='yyyy-mm-dd';
            msg=['  for each path correct all dv/v measurements from the mean dvv taken'];
            msg=[msg,'  between ',datestr(date1,fmt),'-',datestr(date2,fmt)];
            dispc(msg,'c','b');
            I_date=find(h.date2>=date1 & h.date2<= date2);
            ref_mean = mean(h.dvv(I_date,:),'omitnan');
            h.dvv   = h.dvv-ref_mean; 
            npath=numel(ref_mean)   ;
            npath_nan = numel(find(isnan(ref_mean)==1));
            msg=['   => ',num2str(npath_nan),' paths have a nan mean out of ',num2str(npath)];
            msg=[msg,' [',num2str(round(100*npath_nan/npath)),']'];
            dispc(msg,'c','n');
        end
        function add_station_list(h);
            [sta_ Ista]=unique(h.id(:));
            h.sta.lon=double(h.lon(Ista));
            h.sta.lat=double(h.lat(Ista));
            h.sta.sta=h.id(Ista);
        end
        function h=dvv_selection(h_ori,I_path)
            % sort path so that they are by increasing distance:
            [~, J]=sort(h_ori.dvv.dist(I_path));
            I_path=I_path(J);
            % now keep dvv corresponding to I_path:
            h.dvv=h_ori.dvv.dvv(:,I_path);
            h.coh=h_ori.dvv.corrcoef(:,I_path);
            npath=numel(I_path);
            if npath==1
                h.coh=[h.coh h.coh];
                h.dvv=[h.dvv h.dvv];
            end
            % add metadata : 
            h.lon   = h_ori.dvv.lon(:,I_path);
            h.lat   = h_ori.dvv.lat(:,I_path);
            h.id    = h_ori.dvv.id(:,I_path) ; 
            h.dist  = h_ori.dvv.dist(I_path) ;
            h.az    = h_ori.dvv.az(I_path)   ;
            h.baz   = h_ori.dvv.baz(I_path)  ;
            h.in    = h_ori.dvv.in;
            h.date1 = h_ori.dvv.date1 ;
            h.date2 = h_ori.dvv.date2 ;
            h.add_station_list    ;
            h.sta_ori = h_ori.sta ;
            h.fname = h_ori.dvv.file;
            % remove wrong value :
            I_extrem=find(h.dvv(:)==h.in.stretch(1) | h.dvv(:)==h.in.stretch(2) | isnan(h.coh(:))==1);
            h.dvv(I_extrem)=nan;
            h.coh(I_extrem)=nan;
            
            %if exist('titre')
            %    h.titre=titre;
            %end
        end
    end
end