
classdef live < pycorr.live_md & pycorr.live_cvg 
    methods
        function [cc ha] = p_cc_sta_vs_others(h,ksta,p1,p2,varargin)
            in.norm=true   ; 
            in.xlim = []   ; 
            in.flip = false; 
            in.gmap=false  ; 
            %on reporte les options de h.read_cc_sta_vs_others : 
            in.cmp=1     ; % vector of components to be read
            in.day=false ; % vector of date to be read (false = reference)
            in.taper=200 ; % taper length (npts) applied at both side of each cc 
            in.cmap ='gray';
            in.sym = false ;
            in=lang.parse_options(in,varargin); 
            cc=h.read_cc_sta_vs_others(ksta,p1,p2,'cmp',in.cmp,'day',in.day,'taper',in.taper);
            % remove nan lines
            nsta = size(cc.cc,2)
            Inan=[]
            for ista =1 :nsta 
                if numel(   find(isnan(cc.cc(:,ista))==1))>0 | max(cc.cc(:,ista))==0
                    Inan(end+1) = ista
                end
            end
            cc.cc(:,Inan) = []; 
            cc.dist(Inan) = [];
            
            %plot the map of the station :
            ha(1)=subplot(1,2,1);
            in_sta.color=[0.5 0.5 0.5]; 
            if in.gmap 
                h.p_gmap('msize',10,'color','c');
                sta=h.get_station_list; 
                plot(sta.lon(ksta),sta.lat(ksta),'^','color','r','MarkerSize',10,'MarkerFaceColor','r')
	    else
                h.p_map(struct(),struct(),in_sta); 
                h.p_station(ksta,'color','r');
            end
            % plot the correlation vs dist : 
            ha(2) = subplot(1,2,2);
            if in.flip 
                cc.cc = flipud(cc.cc);
            end
            if in.sym 
                cc.cc = cc.cc + flipud(cc.cc);
            end
            if in.norm 
                cc.cc=s2d.norm(cc.cc); 
            end 
            pcolor(h.t,cc.dist,cc.cc');
            if ~isempty(in.xlim) ; 
                xlim(in.xlim);
            end
            colormap(in.cmap)  ;
            shading('interp');
            %get the title of the plot a bit ackward here : 
            aa = str.split(cc.titre{1},'.');
            aa = str.split(aa{end},' ')    ;
            titre={};
            titre(1)={[h.get_station_list.sta{ksta}, ' vs the rest of the network']};
            %if in.norm 
            % 	titre(1)={[titre{1},' (normalized)']};
            % end
            % titre(2)={[aa{end-2},' ',aa{end-1},' ',aa{end}]};
            p.label('correlation time [s]','distance [km]',titre);
            set(gca,'color','w');
        end

        function [st ha]=p_cc_stacked_by_distance(h,p1,p2,varargin)
            in.dist =[0 100] ; 
            in.delta=[2]     ; % +- delta/2
            in.lon  = []     ; % restrict station selection 
            in.lat  = []     ; % to a given area; 
            in.step = 1      ; 
            in.max_path=500  ; %stack pas plus de 500 trajets a chaque fois ?
            in.save = true   ; 
            in.p_npath=true  ;  
            in.xlim=false    ;
            in.sym = false   ;
            in.p_vel=[1,2,3,4]; %add line corresponding to this vel [km/s]
            in.tag='';
            in.cmp=1;
            in= lang.parse_options(in,varargin); 
            %get the number of plot : 
            nplot=1;
            if in.p_npath ; 
                nplot=2;
            end
            %split the figure accordingly : 
            if nplot == 2 
                ha=p.split_matrix_trace_vertical  ;
            else 
                ha =subplot(1,1,1);
            end
            %get the CC stacked : 
            st= get_cc_stacked_by_distance(h,p1,p2,in); 
            % get the xlim and ylim : 
            if ~in.xlim 
                limit= min([h.t(end) round(2*st.dist2(end))]);
                xlimit=[-limit limit];
            else 
                xlimit=in.xlim; 
            end
            ylimit=[st.dist1(1) st.dist(end)];
            
            subplot(ha(1)); 
            pcolor(st.t,st.dist,st.C'); 
            for ivel = 1 :numel(in.p_vel)
                hl=line([0 xlimit(2)],[0 xlimit(2)*in.p_vel(ivel)],'color','r','linewidth',2,'linestyle','--');
            end
            shading('flat') ; 
            xlim(xlimit);
            ylim(ylimit);
            titre = ['cc vs dist ',str.num2str(p1,2),'-',str.num2str(p2,2)];
            titre = [titre,' stacked by distance +-',num2str(st.in.delta),' km'];
            p.label('correlation time [s]','distance [km]',titre,16,'bold',14)
            colormap(gray);
            if in.p_npath 
                subplot(ha(2)); 
                barh(st.dist,st.npath) 
                ylim(ylimit);
                set(gca,'Ytick',[])
                title('number of path stacked','fontweight','bold','fontsize',12)
            end
        end
        
        function cc = p_cc_vs_dist(h,p1,p2,varargin)
            in.norm=true     ; 
            %in.norm_ca=false ; % if norm : should we normalize separately pos/neg time
            in.npath=200  ;
            in.cmp = 1; 
            in.keep_trace_with_data=false;
            in.min_day_of_data=0; 
            in=lang.parse_options(in,varargin)
            
            [dist_ I_dist] = sort(h.dist) ; 
            npath=numel(I_dist);
            J=round(linspace(1,npath,min(npath,in.npath)));
            I_dist=I_dist(J); %indice des correlation a ploter : 
            cc =h.read_cc(I_dist,p1,p2,'cmp',in.cmp);
            
            if in.keep_trace_with_data
                I=find(max(cc.cc>0));
                cc.cc=cc.cc(:,I);
                cc.dist=cc.dist(I);
                cc.ref_nstack = cc.ref_nstack(I);
                cc.az =cc.az(I);
            end
            if in.min_day_of_data >0 
                I=find(cc.ref_nstack>in.min_day_of_data);
                cc.cc=cc.cc(:,I);
                cc.dist=cc.dist(I);
                cc.ref_nstack = cc.ref_nstack(I); 
                cc.az = cc.az(I) ;
            end
            if in.norm 
                cc.cc=s2d.norm(cc.cc);
            end
            pcolor(h.t,cc.dist,cc.cc'); 
            if in.norm
                caxis([-1 1]);
            end
            shading('flat')
            titre=['CC vs dist ',str.num2str(p1,2),'-',str.num2str(p2,2)];
            p.label('correlation time[s]','distance [km]',titre)
        end
                

        function [cc_mat ref ha]=p_matrix(h,kpath,p1,p2,varargin)
        %plot a matrix of daily correlations + the reference
            in.cmp=1    ; 
            in.norm=true;
            in.av=false ; 
            in.lw = 1   ;
            in.xlim=[]  ;
            in=lang.parse_options(in,varargin); 
            
            if isempty(in.xlim)
                xlimit=h.get_xlimit(kpath); 
            else 
                xlimit=in.xlim;
            end
            [cc_mat ref]=h.read_cc_matrix(kpath,p1,p2,in.cmp);
            if in.norm 
                cc_mat=s2d.norm(cc_mat);
            end
            if in.av > 1
                cc_mat=s2d.stack_linear(cc_mat,in.av,1);
            end
            ha=p.split_matrix_trace; 
            %plot matrix : 
            subplot(ha(1))
            imagesc(h.t,mean([h.date1 h.date2]'),cc_mat');
            shading('flat')
            xlim(xlimit)
            datetick('y','yyyy/mm/dd','keepticks','keeplimits')
            ylabel('calendar time','fontweight','bold')
            title(ref.titre,'fontweight','bold');
            %plot the reference : 
            subplot(ha(2))
            plot(h.t,ref.cc,'linewidth',in.lw);
            xlim(xlimit)
            xlabel('correlation time [s]','fontweight','bold','fontsize',14)
            grid minor
        end
        

        function [cc ha]=p_cc1(h,kpath,p1,p2,varargin)
        %plot a single correlation with map 
            ha(1)=subplot('Position',[0.1 0.5 0.8 0.45]);
            h.p_station_pair(kpath);
            ha(2)=subplot('Position',[0.1 0.1 0.8 0.25]) ;
            if numel(varargin) ~=0 
                arg={varargin{1,:},'ha',ha(2)};
            else
                arg={'ha',ha(2)};
            end
            [cc]=h.p_cc(kpath,p1,p2,arg{1,:});
        end
	  

        function [cc ha]=p_cc(h,kpath,p1,p2,varargin)
        %generic function to plot a set of correlations. 
            in.day=false ;  
            in.cmp=1     ;  
            in.taper=200 ;	
            in.norm=true ;
            in.snr=true  ;
            in.color='k' ;
            in.lw=1      ;     
            in.dy=0      ;
            in.ttl_fs =12;     
            in.ttl_leg=false ; 
            in.collapse=false;
            in.grid='on'  ;
            in.xlim = []  ;
            in.ha=false   ; %used to force plot in an axes (ok only if single trace or collapse)
            in=lang.parse_options(in,varargin);
            
            %read all correlations  :
            cc =h.read_cc(kpath,p1,p2,'day',in.day,'cmp',in.cmp);
            ncc = size(cc.cc,2);
            %signal processing : 
            if in.norm 
                cc.cc=s2d.norm(cc.cc);
            end
            %determine xlimit :
            xlimit=h.get_xlimit(kpath);
            if ~isempty(in.xlim)
                xlimit=in.xlim ;
            end
            %determine the color of the plot : 
            if ncc > 1 & numel(in.color)==1 
                in.color=char(double(in.color*ones(1,ncc)));
            end
            %determine the linewidth :
            if ncc > 1 & numel(in.lw)==1
                in.lw= in.lw*ones(1,ncc);
            end
            
            %determine how many plot do we have and create subplots: 
            if in.ha == false 
                if in.collapse==true 
                    xx=p.split_horizontal(1,'dy',0) ; 
                    ha=ones(1,ncc);
                    ha(1:end)=xx ;
                    hold on
                else
                    ha=p.split_horizontal(ncc,'dy',in.dy) ; 
                end
            else 
                hold on
                ha=ones(1,ncc);
                ha(1:end)=in.ha;
            end
            %lets plot : 
            for icc =1 : ncc 
                subplot(ha(icc))
                plot(h.t,cc.cc(:,icc),'color',in.color(icc),'linewidth',in.lw(icc));
                if in.snr 
                    snr=cc2d.sw_snr(cc.cc(:,icc),[],h.t,cc.dist(icc),struct('p1',cc.p1(icc),'p2',cc.p2(icc)));
                    snr_a=str.num2str(log(snr.amp(2)/snr.noiz(2)),4);
                    snr_c=str.num2str(log(snr.amp(1)/snr.noiz(1)),4);
                    cc.titre(icc)={[cc.titre{icc},'  SNR ',snr_a,'-',snr_c]};
                    hold on
                    p_rsb(snr,h.t);
                end	
                set(gca,'XminorGrid',in.grid);
                set(gca,'YminorGrid',in.grid);
                %end
                xlabel('');
                xlim(xlimit);
            end		
            xlabel('correlation time [s]','fontweight','bold','fontsize',12);
            
            %add title/legend : 
            if in.collapse==true
                hl=legend(cc.titre);
                legend('boxoff') 
                set(hl,'fontweight','bold','fontsize',in.ttl_fs);
            else 
                for icc=1:ncc
                    subplot(ha(icc))
                    if in.ttl_leg==false 
                        title(cc.titre{icc},'fontweight','bold','fontsize',in.ttl_fs)
                    else 
                        hl=legend(cc.titre{icc},'fontweight','bold','fontsize',in.ttl_fs);
                        legend('boxoff') ;
                        set(hl,'fontweight','bold','fontsize',in.ttl_fs);
                    end	
                end
            end
        end
        
        %-------------------------------------------------------------------
        %                      UTILITY FUNCTIONS 
        %------------------------------------------------------------------
        function bs= get_cc_stacked_by_distance(h,p1,p2,varargin)
            in.dist =[0 100] ; 
            in.delta=[2]     ; % +- delta/2
            in.lat  =[]      ; 
            in.lon  =[]      ;
            in.step = 1      ; 
            in.max_path=500  ; % stack pas plus de 500 trajets a chaque fois ?
            in.save = false  ; % save the result in a file ?
            in.tag= '';
            in.cmp=1;
            in.max_pos=true  ; %put the maximum in the positive time ?
            in= lang.parse_options(in,varargin); 
            limit_coord = true    ;
            if isempty(in.lat)
                in.lat = [-90 90]   ;
                limit_coord = false ;
            end
            if isempty(in.lon) 
                in.lon=[-360 360]   ; 
                limit_coord = false ;
            end
            %determine distance interval on which we will stack: 
            in.dist(2) = min([in.dist(2)  max(h.dist)]); 
            dist = [in.dist(1):in.step:in.dist(2)];
            dist1 = dist- in.delta/2; 
            dist2 = dist+ in.delta/2; 
            dist1(dist1<0)=0;
            ndist = numel(dist1); 
            
            stack = zeros(1,ndist) ; 
            %get indice of the distance vector where we have indeed some data : 
            for idist = 1 : ndist 
                stack(idist) = numel(find(h.dist >= dist1(idist) & h.dist <= dist2(idist)));
            end
            I_data = find(stack > 0) ; 
            %get the data and the stack : 
            ndist_data = numel(I_data)  ; 
            nwf = numel(h.t)            ; 
            bs = struct()               ;
            bs.p1=p1                    ;
            bs.p2=p2                    ;			
            bs.C = zeros(nwf,ndist_data);
            bs.dist1 = dist1(I_data)    ;
            bs.dist2 = dist2(I_data)    ; 
            bs.dist  = dist(I_data)     ; 
            bs.npath = zeros(1,ndist_data);
            bs.t     = h.t; 
            bs.tau   = h.tau;
            bs.in    = in ; 
            bs.dir   = 'stack_by_dist'; 
            bs.fname = [bs.dir,'/st_',in.tag,'_',str.num2str(p1,5),'-',str.num2str(p2,5),'s'];
            if limit_coord 
                bs.fname =[bs.fname,'_coord_',num2str(in.lon(1)),'_',num2str(in.lon(2))] ;
                bs.fname =[bs.fname,'_',num2str(in.lat(1)),'_',num2str(in.lat(2))] ; 
            end
            bs.fname = [bs.fname,'__dist_',num2str(round(in.dist(1))),'_',num2str(round(in.dist(2))),'km'];
            bs.fname = [bs.fname,'__step_',num2str(in.step),'__delta_',num2str(in.delta)];
            if in.sym 
                bs.fname =[bs.fname,'_sym'];
            end
            bs.fname = [bs.fname,'__max_path_',num2str(in.max_path),'.mat'];
            
            % 1st get indice of the path having the right coordinante : 
            if numel(in.lat) == 2 % rectangle x/y 
                I_lat1 = find(h.lat(1,:) >= in.lat(1) & h.lat(2,:) >= in.lat(1));
                I_lat2 = find(h.lat(1,:) <= in.lat(2) & h.lat(2,:) <= in.lat(2));
                I_lon1 = find(h.lon(1,:) >= in.lon(1) & h.lon(2,:) >= in.lon(1));
                I_lon2 = find(h.lon(1,:) <= in.lon(2) & h.lon(2,:) <= in.lon(2));
                I_coor = intersect(I_lat1, I_lat2) ; 
                I_coor = intersect(I_coor,I_lon1)  ; 
                I_coor = intersect(I_coor,I_lon2)  ; 
            else 
                I_coor=inpolygon(h.lon,h.lat,in.lon,in.lat);
                I_coor = find(I_coor(1,:) == 1 & I_coor(2,:)==1);
            end
            bs.I_coor = I_coor; 
            % loop on each distance intervall : 
            for idist = 1 : ndist_data	
                I_path = find(h.dist >= bs.dist1(idist) & h.dist <=bs.dist2(idist)); 
                I_path = intersect(I_path,I_coor); 
                if numel(I_path) > in.max_path 
                    I_path = I_path(1:in.max_path); 
                end
                bs.npath(idist)=numel(I_path); 
                % read the all correlation, and put the part with best snr in the pos time:
                c1=h.read_cc(I_path,p1,p2,'cmp',in.cmp);
                if in.max_pos == true 
                    for icc = 1 :numel(c1.dist); 
                        snr=cc2d.sw_snr(c1.cc(:,icc),[],h.t,c1.dist(icc),struct('p1',c1.p1(icc),'p2',c1.p2(icc)));
                        snr_a=log(snr.amp(2)) ;%/snr.noiz(2));
                        snr_c=log(snr.amp(1)) ;%/snr.noiz(1));
                        if snr_a > snr_c 
                            c1.cc(:,icc) = flipud(c1.cc(:,icc));
                        end
                    end
                end
                % now stack and store for this dist interval : 
                bs.C(:,idist)=s2d.norm(sum(s2d.norm(c1.cc),2));
            end
            if in.sym 
                bs.C = bs.C + flipud(bs.C);
            end
            %remove nan rows :
            ndist = size(bs.C,2);
            Inan = [] ;
            for idist = 1 : ndist 
                if numel(find(isnan(bs.C(:,idist))==1)) > 0 | max(bs.C(:,idist)) ==0;
                    Inan(end+1)= idist;
                end
            end
            bs.C(:,Inan)=[];
            bs.dist1(:,Inan)=[];
            bs.dist2(:,Inan)=[];
            bs.dist(:,Inan)=[];
            bs.npath(:,Inan)=[];
            if in.save 
                os.mkdir(bs.dir); 
                save(bs.fname,'-struct','bs');
            end
        end
        
        function xlimit=get_xlimit(h,kpath)
            dist_=max(h.dist(kpath));	
            if dist_ > 100 ;
                limit= min([h.t(end) round(2*dist_)]);
                xlimit=[-limit limit];
            else
                limit =min([h.t(end) 200]);
                xlimit=[-limit limit]; 
            end	
        end	
        x
        function h=live(in)
            h = h@pycorr.live_md(in);
        end
    end
end



%-----------------------------------------------------------------
%                         UTILITY FUNCTIONS USED TO PLOT 
%-----------------------------------------------------------------

function p_rsb(snr,t)
    in.dy=0;
    for ica=1:2
  	xbox=t([snr.I1(ica) snr.I2(ica) snr.I2(ica) snr.I1(ica)]);
        ybox=[-1+in.dy -1+in.dy -0.7+in.dy -0.7+in.dy];
        patch(xbox,ybox,'r','faceAlpha',0.3,'EdgeColor','r');
        xbox=t([snr.Isurf1(ica) snr.Isurf2(ica) snr.Isurf2(ica) snr.Isurf1(ica)]);
        ybox=[0.7+in.dy 0.7+in.dy -0.7+in.dy -0.7+in.dy];
        patch(xbox,ybox,'b','faceAlpha',0.2,'EdgeColor','b');
        xbox=t([snr.Inoise1(ica) snr.Inoise2(ica) snr.Inoise2(ica) snr.Inoise1(ica)]);
        patch(xbox,ybox,'y','faceAlpha',0.2,'EdgeColor','y');
        %if ica==1 & in.decoration
        %	legend('looking for SW','SW identified','noise');
        %end
    end
end	


