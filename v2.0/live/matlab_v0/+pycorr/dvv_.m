%TODO : ameliorer calcul rms :) 
%reecrire patch
%
% similar path : h.dvv.id{10} <=>  h.id(:,h.dvv.I(10))
%
% Cette classe permet de manipuler l'ensemble des mesures de dv/v présent dans un ficher de résultat.
%   - lire les données 
%   - plot dvv moyen 
%   - carte des station, plot des correlation via pycorr.live => ici on est encore en lien avec les CC 
%   - selectionner un sous-ensemble de trajet, et appelle dvv_selection pour manipuler ce sous-ensemble
%     (plot/carte des dv/v, re-selection des dv/v ayant les bon cc, les trajets ayant suffisament de
%     donnees etc, recalage des courbes,...) 
%

classdef dvv_ < pycorr.live
    properties
        dvv=struct()
        sta=[];
        I_id=[];
    end
    methods
        function out=dvv_build_video_around_station(h,varargin)
            in.dist_area=[ 50 100 150 200 250 300 400 500]  ;
            in.dist =[0 inf] ; % distance interstation minimum/maximum
            in.out_dir='Cvideo';
            in.min_day_per_path=400;
            in.cc=[0.5 1];
            in.min_path=30;
            in.ref_period=[]; %h.dvv.date2(1) h.dvv.date2(end)];
            in = lang.parse_options(in,varargin);
            
            os.mkdir(in.out_dir);
            %initialize output struct :
            ndate = numel(h.dvv.date2)     ;
            nsta = numel(h.sta.sta)        ;
            out=struct()                   ;
            out.dvv      = nan(ndate,nsta) ;
            out.coh      = nan(ndate,nsta) ;
            out.dvv_std  = nan(ndate,nsta) ;
            out.rms      = nan(ndate,nsta) ;
            out.date2    = h.dvv.date2     ; 
            out.date1    = h.dvv.date1     ;
            out.dist_area= nan(1,nsta)     ;
            out.npath    = nan(1,nsta)     ; 
            out.in       = in              ;
            out.in_dvv   = h.dvv.in        ;
            out.in_cc    = h.in            ;
            out.sta      = h.sta           ;
            % determine output file name : 
            aa=strsplit(h.dvv.file,'/');
            out.name =[in.out_dir,'/',aa{2}(1:end-3)];
            out.name =[out.name,'__dist_area1_',num2str(in.dist_area(1)),'km'];
            out.name =[out.name,'__min_day_per_path_',num2str(in.min_day_per_path)];
            out.name =[out.name,'__min_path_',num2str(in.min_path)];
            out.name =[out.name,'.mat'];
            if exist(out.name,'file')
                dispc(['  ',out.name,' already exist'],'r','b')
                return
            end
            for ista = 1 :nsta
                for dist_area=in.dist_area
                    [h2 titre]=h.dvv_p_around_station(ista,'dist_area',dist_area,'dist',in.dist,'plot',false);
                    if size(h2.dvv,2) > in.min_path
                        break
                    end
                end
                h2.rm_path_with_less_than_nday_of_data(in.min_day_per_path) ;
                if numel(in.ref_period>0)
                    h2.correct_each_dvv_from_mean(in.ref_period(1),in.ref_period(2)); %date reseau temporaire
                end
                h2.keep_measurement_with_cc_btw(in.cc)                          ;
                out.dvv(:,ista)=mean(h2.dvv,2,'omitnan');
                out.coh(:,ista)=mean(h2.coh,2,'omitnan');
                out.npath(ista)=size(h2.dvv,2);
                out.dist_area(ista)=dist_area;
                out.dvv_std= std(h2.dvv,[],2,'omitnan');
                out.rms(:,ista)=h2.compute_rms ; % ok tant qu'on a qu'une seule cmp
            end
            save(out.name,'-struct','out')
        end
        
        function [h2 titre]=dvv_p_around_station(h,ksta,varargin)
            in.dist_area=[1000]   ; %zone autour de la station
            in.dist =[0 inf]      ; % distance interstation minimum/maximum
            in.plot=true          ;
            in= lang.parse_options(in,varargin);

            lon = h.sta.lon(ksta);
            lat = h.sta.lat(ksta); 
            lat_d = in.dist_area/m_idist(lon,lat,lon,lat+1)*1000;
            lon_d = in.dist_area/m_idist(lon,lat,lon+1,lat)*1000;
            
            % la station 1 est dans la bonne gamme de longitude/latitude :
            I1=find(h.dvv.lon(1,:) >=lon-lon_d & h.dvv.lon(1,:)<= lon+lon_d);
            I2=find(h.dvv.lat(1,:) >=lat-lat_d & h.dvv.lat(1,:)<= lat+lat_d);
            % la station 2 est dans la bonne gamme de longitude/latitude :
            I3=find(h.dvv.lon(2,:) >=lon-lon_d & h.dvv.lon(2,:)<= lon+lon_d);
            I4=find(h.dvv.lat(2,:) >=lat-lat_d & h.dvv.lat(2,:)<= lat+lat_d);
            % les deux stations sont dans les bonnes zonnes : 
            I_path=intersect(I4,intersect(I3,intersect(I1,I2)));
            % on ne garde que les trajets ayant les bonnes distances : 
            J_path=find(h.dvv.dist(I_path)>= in.dist(1) & h.dvv.dist(I_path)<=in.dist(2));
            I_path=I_path(J_path);

            titre2 = ['Area of ',num2str(in.dist_area),'km'];
            titre2 = [titre2,'  inter-station dist ',num2str(in.dist(1)),'-',num2str(in.dist(2)),'km'];
            titre = {[h.sta.sta{ksta},' : ', h.dvv_get_title_generic]};
            titre(2)={titre2};
            
            h2 = pycorr.dvv_selection(h,I_path);
            
            %titre2 = ['Area of ',num2str(in.dist_area),'km'];
            %titre2 = [titre2,'  inter-station dist ',num2str(in.dist(1)),'-',num2str(in.dist(2)),'km'];
            %titre = {[h.sta.sta{ksta},' : ', h.dvv_get_title_generic]};
            %titre(2)={titre2};
            if in.plot
                ha=h2.plot_map_and_matrix;
                subplot(ha(1))
                title(titre,'fontweight','bold','fontsize',14);
            end
            return
            %out=h.read_and_clean_matrix_for_selected_path(I_path,in);
            %ha=h.dvv_out_p_map_and_matrix(out);
            %add title
            %subplot(ha(1))
            %titre2 = ['Area of ',num2str(in.dist_area),'km'];
            %titre2 = [titre2,'  inter-station dist ',num2str(in.dist(1)),'-',num2str(in.dist(2)),'km'];
            %titre = {[h.sta.sta{ksta},' : ', h.dvv_get_title_generic]};
            %titre(2)={titre2};
            %title(titre,'fontweight','bold','fontsize',14);
        end
        function [h2]=dvv_p_single_station(h,ksta,varargin)
            in.dist = [0 inf] ;
            in.xlimit=false;
            in.plot=false;
            in = lang.parse_options(in,varargin); 
            
            %first find measurements around the station ksta: 
            I=find(h.dvv.I_id(1,:)==ksta | h.dvv.I_id(2,:) ==ksta);
            J=find(h.dvv.dist(I)>=in.dist(1) & h.dvv.dist(I) <=in.dist(2));
            I_path=I(J);
            h2 = pycorr.dvv_selection(h,I_path);
        end
        function h2 = dvv_select_path_in_area(h,lon,lat,varargin)
            in = struct() ; 
            in.dist =[0 inf] ; 
            in = lang.parse_options(in,varargin);  
            
            I_lon1 = find(h.dvv.lon(1,:) >= lon(1) & h.dvv.lon(2,:)>=lon(1));
            I_lon2 = find(h.dvv.lon(1,:) <= lon(2) & h.dvv.lon(2,:)<=lon(2));
            I_lat1 = find(h.dvv.lat(1,:) >= lat(1) & h.dvv.lat(2,:)>=lat(1));
            I_lat2 = find(h.dvv.lat(1,:) <= lat(2) & h.dvv.lat(2,:)<=lat(2));
            
            I_path=intersect(I_lon1,I_lon2);
            I_path=intersect(I_path,I_lat1);
            I_path=intersect(I_path,I_lat2);
            
            I_dist = find(h.dist >= in.dist(1) & h.dist <=  in.dist(2)); 
            I_path = intersect(I_dist,I_path);
            
            h2 = pycorr.dvv_selection(h,I_path);
        end
        function correct_each_dvv_from_mean(h,date1,date2)
            fmt='yyyy-mm-dd';
            msg=['  for each path correct all dv/v measurements from the mean dvv taken'];
            msg=[msg,'  between ',datestr(date1,fmt),'-',datestr(date2,fmt)];
            dispc(msg,'c','b');
            I_date=find(h.dvv.date2>=date1 & h.dvv.date2<= date2);
            ref_mean = mean(h.dvv.dvv(I_date,:),'omitnan');
            h.dvv.dvv   = h.dvv.dvv-ref_mean; 
            npath=numel(ref_mean)   ;
            npath_nan = numel(find(isnan(ref_mean)==1));
            msg=['   => ',num2str(npath_nan),' paths have a nan mean out of ',num2str(npath)];
            msg=[msg,' [',num2str(round(100*npath_nan/npath)),']'];
            dispc(msg,'c','n');
        end

        function keep_measurement_with_cc_btw(h,cc)
            I=find(h.dvv.corrcoef(:)< cc(1) | h.dvv.corrcoef(:) > cc(2));
            h.dvv.dvv(I)=nan;
            h.dvv.corrcoef(I)=nan;
            cc_str=[num2str(cc(1)),'-',num2str(cc(2))];
            msg1=['  setting all dv/v measurements with a cc outside ',cc_str,' to nan'];
            msg2=['  => ',num2str(numel(I)),' meas out of ',num2str(numel(h.dvv.dvv(:))),' set to nan'];
            dispc(msg1,'c','b');
            dispc(msg2,'c','n');
        end
        function rm_path_with_less_than_nday_of_data(h,nday)
          % remove path that have less than nday of measurements (we correct the numbef of meas
            % by stackshift so that it is in day
            nday_per_path = numel(h.dvv.date2)-sum(isnan(h.dvv.corrcoef));
            nday_per_path = nday_per_path*double(h.dvv.in.stack_shift);
            I_path=find(nday_per_path<nday);
            h.rm_path(I_path);
            msg1=['  removing all path with less than ',num2str(nday),' days of data'];
            msg2=['  => ',num2str(numel(I_path)),' paths removed out of ',num2str(numel(nday_per_path))];
            dispc(msg1,'c','b');
            dispc(msg2,'c','n');
        end
        function rm_path(h,I_path)
            h.dvv.dvv(:,I_path)=[];
            h.dvv.corrcoef(:,I_path)=[];
            h.dvv.lon(:,I_path)=[];
            h.dvv.lat(:,I_path)=[];
            h.dvv.id(:,I_path)=[];
            h.dvv.dist(I_path)=[];
            h.dvv.az(I_path)=[];
            h.dvv.baz(I_path)=[];
            h.dvv.I_id(:,I_path)=[];
            h.dvv.I(I_path)=[];
        end
        function titre=dvv_get_title_generic(h)
            titre=[];
            titre = ['dv/v ',num2str(h.dvv.in.p1),'-',num2str(h.dvv.in.p2),'s'];
            titre = [titre,'   stack ',num2str(h.dvv.in.stack),'-',num2str(h.dvv.in.stack_shift)];
            titre = [titre,'   coda length ',[num2str(h.dvv.in.tf),'s']];
        end
        function [x ha]=dvv_p_mean(h,varargin)
            in=struct()
            in.dist=[0 inf];
            in.cc_min=0    ;
            if isfield(h.dvv.in_ml,'cmp')
                in.cmp =[1: numel(h.dvv.in_ml.cmp)] ; 
            else 
                in.cmp='Z';
            end
            in.p_rms=false  ;
            in.p_nmeas=false;
            in.p_mcc  =false; 
            in.color='b'    ;
            in.herror=false ; 
            in=lang.parse_options(in,varargin);
            
            nplot=1+in.p_nmeas + in.p_mcc;
            if nplot >1
                ha=p.split_horizontal(nplot,'dy',0.0);
                subplot(ha(1));
            end
            hold on
            x=h.dvv_select_path('dist',in.dist,'cc_min',in.cc_min,'cmp',in.cmp);
            if in.herror == false 
                if numel(x.cmp_dvv) > 1 % depends si il y'a une ou plusieurs composante 
                    plot(h.dvv.date,x.cmp_dvv,'-','color',in.color);
                else 
                    plot(h.dvv.date,x.dvv_mean,'-','color',in.color);
                end
            else
                herror=(h.dvv.date2-h.dvv.date1)/2;
                p.herrorbar(h.dvv.date,x.cmp_dvv,herror,herror,'-ob')%,'-*','color',in.color);
                p.herrorbar(h.dvv.date,x.cmp_dvv,herror,herror,'.r')%,'*','color','r');
            end
            datetick('x','keeplimits');
            if nplot==1
                xlabel('calendar time','fontweight','bold','fontSize',14);
            end
            ylabel('dv/v','fontweight','bold','fontSize',14);
            if in.p_rms
                a=[h.dvv.date h.dvv.date(end:-1:1)];
                b=[x.cmp_dvv+x.cmp_rms ; x.cmp_dvv(end:-1:1)-x.cmp_rms(end:-1:1)];
                patch([h.dvv.date h.dvv.date(end:-1:1)],[x.cmp_dvv+x.cmp_rms ; ...
                                    x.cmp_dvv(end:-1:1)-x.cmp_rms(end:-1:1)]', ...
                      in.color,'facealpha',0.2,'EdgeColor',[1 1 1]);
            end
            %..
            last_subplot=1;
            %..
            if in.p_mcc
                last_subplot=last_subplot+1;
                subplot(ha(last_subplot))
                plot(h.dvv.date,x.cc_mean);
                datetick('x','keeplimits');
                ylabel('mean cc','fontweight','bold','fontSize',14);
            end
            %..
            if in.p_nmeas
                last_subplot=last_subplot+1;
                subplot(ha(last_subplot));
                if numel(x.cmp_nmeas) > 1
                    bar(h.dvv.date,x.cmp_nmeas,'b');
                else
                    bar(h.dvv.date,x.nmeas,'b')
                end
                datetick('x','keeplimits');
                ylabel('number of paths','fontweight','bold','fontSize',14);
            end
            xlabel('calendar time','fontweight','bold','fontSize',14);
        end
        
        
        function x=dvv_select_path(h,varargin)
        % select paths having the appropriate distance/cc_min 
            in.dist=[20 100]                    ;
            in.cc_min=0.8                       ;
            if isfield(h.dvv.in_ml,'cmp')
                in.cmp =[1: numel(h.dvv.in_ml.cmp)] ; 
            else 
                in.cmp='Z';
            end
            in=lang.parse_options(in,varargin)  ;
            
            [ndate npath] =size(h.dvv.dvv);
            ncmp = numel(in.cmp)          ;
            I_dist=find(h.dist<=in.dist(2) & h.dist >= in.dist(1) );
            x.dvv=nan(ndate,npath,ncmp) ;
            x.dvv_mean=zeros(ndate,ncmp);
            x.dvv_med =zeros(ndate,ncmp);
            x.cc_mean =zeros(ndate,ncmp);
            x.nmeas   =zeros(ndate,ncmp);
            
            for icmp=1:ncmp 
                for idate=1:ndate
                    I_cc=find(h.dvv.corrcoef(idate,:,icmp) >= in.cc_min);
                    I=intersect(I_dist,I_cc);
                    if ~isempty(I)
                        x.dvv_mean(idate,icmp)= mean(h.dvv.dvv(idate,I,icmp))     ;
                        x.dvv_med(idate,icmp) = median(h.dvv.dvv(idate,I,icmp))   ;
                        x.cc_mean(idate,icmp) = mean(h.dvv.corrcoef(idate,I,icmp));
                        x.nmeas(idate,icmp)   = numel(I)                          ;
                        x.dvv(idate,I,icmp)   = h.dvv.dvv(idate,I,icmp)           ;
                    end
                end
                x.cmp_dvv  =mean(x.dvv_mean');
                x.cmp_cc   =mean(x.cc_mean') ;
                x.cmp_nmeas=mean(x.nmeas')   ;
            end
            %now : compute rms value : 
            x.rms= zeros(ndate,ncmp)             ;
            wc=1./mean([h.dvv.in.p1 h.dvv.in.p2]); %Hz
            T= 1./(h.dvv.in.p2-h.dvv.in.p1)      ; %T [Hz]
            dist_mean=mean(h.dist(I_dist))       ;
            t1_mean=dist_mean/h.dvv.in.v1+h.dvv.in.coda_dp*h.dvv.in.p2;
            t2_mean=t1_mean+h.dvv.in.tf;%h.mt.it(1).tf;
            for icmp =1 : ncmp 
                rms1= sqrt((6*T*sqrt(pi/2))/(wc^2*(t2_mean^3-t1_mean^3)))  ;
                rms2= (sqrt(1-x.cc_mean(:,icmp).^2))./(2*x.cc_mean(:,icmp));
                x.rms(:,icmp)=rms1*rms2./sqrt(x.nmeas(:,icmp))             ;
            end
            x.cmp_rms=mean(x.rms');
        end
        function add_I_id(h)
          % pour accelerer la recherche d'indice, ajoute une matrice d'indice de station 
          % [2 ncpl] : 
            ncpl = size(h.dvv.id,2);
            nsta = numel(h.sta.sta);
            h.dvv.I_id = zeros(2,ncpl);
            dispc('  here we use only 10 character to build h.dvv.I_id->to be removed in the futur','g','n')
            for ista = 1 : nsta
                I=find(strncmp(h.dvv.id(1,:),h.sta.sta{ista},10)==1);
                h.dvv.I_id(1,I)=ista;
                I=find(strncmp(h.dvv.id(2,:),h.sta.sta{ista},10)==1);
                h.dvv.I_id(2,I)=ista;
            end
            %--
            %tic
            %ncpl = size(h.id,2);
            %I_id2 = zeros(2,ncpl);
            %nsta = numel(h.sta.sta);
            %for icpl =1 : ncpl 
            %    I=find(strcmp(h.id(1,icpl),h.sta.sta));
            %    I_id2(1,icpl)=I; 
            %    I=find(strcmp(h.id(2,icpl),h.sta.sta));
            %    I_id2(2,icpl)=I; 
            %end
        end
        
        function h=dvv_(pydb_file)
            h = h@pycorr.live(pydb_file)                                       ;
            h.dvv=h5.read_group(pydb_file,'/dvv');
            try
                [h.dvv.coeffmatrix]        = h5.tree_to_matrix(pydb_file,'/dvv/coeffmatrix');
            catch
                a=1;
            end
            [h.dvv.dvv dset_path]      = h5.tree_to_matrix(pydb_file,'/dvv/dvv')        ;
            [h.dvv.corrcoef]           = h5.tree_to_matrix(pydb_file,'/dvv/corrcoef')   ;

            h.dvv.in          = h5.read_group(pydb_file,'/dvv/in_')            ;
            h.dvv.in_ml       = h5.read_group(pydb_file,'/dvv/in_ml')          ;
            h.dvv.md_c        = h5.read_group(pydb_file,'/dvv/md_c')           ;
            h.dvv.date1 = h.dvv.md_c.date1+366 ;
            h.dvv.date2 = h.dvv.md_c.date2+366 ;
            h.dvv.date  = h.dvv.md_c.date +366 ;
            h.dvv.in.tf = double(h.dvv.in.tf);
            h.sta = h.get_station_list();
            % mapping orrelation metadata and dvv  :
            h_id=cell_.merge(h.id(1,:),h.id(2,:),'/');
            npath = numel(dset_path);
            h.dvv.I =zeros(1,npath);
            h.dvv.id = cell(2,npath);
            %h.dvv.id = dset_path;
            %..
            dispc(' to overcome the id bug, compare only id on 21 caractere','g','n')
            for ipath = 1 : numel(dset_path)
                dispc([str.num2str(ipath,5),' : ' , str.num2str(numel(dset_path{ipath}),4)],'b','b');
                aa=strsplit(dset_path{ipath},'/');
                h.dvv.id(1,ipath)=aa(1);
                h.dvv.id(2,ipath)=aa(2);
                try 
                    h.dvv.I(ipath)=find(strncmp(dset_path{ipath},h_id,21)==1);
                    %h.dvv.I(ipath)=strmatch(dset_path{ipath},h_id);
                catch	
                    keyboard
                end
            end	
            h.dvv.dist = h.dist(h.dvv.I);
            h.dvv.az   = h.az(h.dvv.I)  ; 
            h.dvv.baz  = h.baz(h.dvv.I) ;
            h.dvv.lon  = h.lon(:,h.dvv.I);
            h.dvv.lat  = h.lat(:,h.dvv.I);
            h.add_I_id;
            h.dvv.file = pydb_file;
        end
    end  
end


