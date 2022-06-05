% Generate a db_disp*file from the h5 containing only the metadata to speedup everything 
% if ~pymd_disp*mat : load pydb and generate the pymd_disp*mat
% - stream matrix & disp directly from disk. 
% - live_md : read  CC from pydb_disp*h5 file or db_disp*mat file 
% on numerote les trajets comme d'apres l'ordre des CC pas des courbes de dispersion !! 
% => kpath = 1 => 1ere CC pas la premiere courbe de dispersion !! 
% h.disp.pydb_file : is always the name of the .h5 files containing everything ! 
%
% acausal, causal, sym 
classdef disp_ < pycorr.live 
    properties
        disp=struct(); 
    end
    methods 
        function disp_concatenate_meas_at_period(h,kper,kcmp,varargin)
            in.sym = false ; 
            in = lang.parse_options(in,varargin);
            % for a single period : summarize all measurements : SNR, disp c/a, dist, az, baz, id .... 
            % save in a single mat file for further processing ? 
            per      = h.disp.per(kper)    ; 
            aa= str.split(h.disp.pydb_file);
            out_dir = [aa{1},'/tomo_',aa{2}(5:end-3)]              ;
            out_file= [out_dir,'/meas_',str.num2str(per,3),'s.mat'];
            os.mkdir(out_dir)        ;
            npath = numel(h.disp.dist);
            E = struct()          ; 
            E.id  = cell(2,npath) ;
            E.dist= zeros(1,npath); 
            E.az  = zeros(1,npath);
            E.baz = zeros(1,npath);
            E.lon = zeros(2,npath); 
            E.lat = zeros(2,npath); 
            E.gv  = zeros(2,npath); 
            E.snr = zeros(2,npath); 
            E.per = per           ; 
            E.in_file = h.disp.pydb_file;
            E.file    = out_file        ; 
            for ipath = 1 : npath
                disp_= h.disp_read(ipath,kcmp); 
                E.id(:,ipath)  = disp_.id ; 
                E.lat(:,ipath) = disp_.lat;
                E.lon(:,ipath) = disp_.lon;
                E.dist(ipath) = disp_.dist;
                E.az(ipath)   = disp_.az  ;
                E.baz(ipath)  = disp_.baz ;
                E.gv(:,ipath) = disp_.gv(kper,1:2) ; 
                E.snr(:,ipath)= disp_.snr(kper,1:2);
            end
            save(out_file,'-struct','E');
        end
    end
    methods % plotting methods 
        function disp_p_matrix(h,kpath,kcmp,varargin)
            in.p1 = 5;  % filter used to display the cc 
            in.p2 = 20; 
            in.caxis=[0 10];
            in.p_trace = false ;
            in.p1 = 5;  % to filte the trace 
            in.p2 =20;  % :-) 
            in = lang.parse_options(in,varargin);
            % split the figure according to users wish :)
            width=0.4;
            if in.p_trace 
                ha(1)= subplot('position',[0.1 0.4 width 0.5]);
                ha(2)= subplot('position',[0.1+width+0.05 0.4 width 0.5]);
                ha(3)= subplot('position',[0.1 0.1 0.7750 0.20]) ;
            else 
                ha(1)= subplot('position',[0.1 0.11 width 0.815]);
                ha(2)= subplot('position',[0.1+width+0.05 0.11 width 0.815]);
            end
            % read the data 
            disp_ = h.disp_read(kpath,kcmp,true) ; 
            % plot acausal part : 
            subplot(ha(1))
            h.disp_pa_matrix(disp_,1,'pos',get(ha(1),'Position'),'colorbar',false);			
            % plot causal : 
            subplot(ha(2))
            h.disp_pa_matrix(disp_,2,'pos',get(ha(2),'Position'),'colorbar',false);
            % plot cc waveform : 
            if in.p_trace 
                h.p_cc(h.disp.I(kpath),in.p1,in.p2,'ha',ha(3));
                title('')
            end
            % add title : 
            subplot(ha(1));
            titre = h.disp_get_title(kpath,kcmp) ; 
            title(titre,'fontweight','bold');
        end
        
        function disp_p(h,kpath,kcmp,varargin)
            in.p1 = 5;  % filter used to display the cc 
            in.p2 = 20; 
            in.caxis=[0 10];
            in.p_trace = false ;
            in = lang.parse_options(in,varargin);
            % split the figure according to users wish :)
            width=0.37;
            if in.p_trace 
                ha(1)= subplot('position',[0.1 0.4 width 0.5]);
                ha(2)= subplot('position',[0.1+width+0.05 0.4 width 0.5]);
                ha(3)= subplot('position',[0.1 0.1 0.7750 0.20]) ;
            else 
                ha(1)= subplot('position',[0.1 0.11 width 0.815]);
                ha(2)= subplot('position',[0.1+width+0.05 0.11 width 0.815]);
            end
            % read the data 
            disp_ = h.disp_read(kpath,kcmp,false) ; 
            % plot acausal part : 
            subplot(ha(1))
            h.disp_pa_disp(disp_,1);
            grid minor 
            colormap(jet)
            caxis(in.caxis);
            % plot causal : 
            subplot(ha(2))
            h.disp_pa_disp(disp_,2)
            caxis(in.caxis)
            grid minor 
            hc=colorbar; 
            aa=get(hc,'Position');
            set(hc,'Position',[0.91 aa(2) 0.01 aa(4)])
            % plot cc waveform : 
            if in.p_trace 
                h.p_cc(h.disp.I(kpath),in.p1,in.p2,'ha',ha(3));
                title('')
            end
            % add title : 
            subplot(ha(1));
            titre = h.disp_get_title(kpath,kcmp) ; 
            title(titre,'fontweight','bold');
        end
    end
    methods % atomic plotting method : do not init axex and dont read data.
        function disp_pa_matrix(h,disp_,kca,varargin)
        % plot the energy diagram + dispersion curve 
            in.color='k'  ;  % disp curve color 
            in.lw   = 1   ;  % linewidth of the disp curve 
            in.s    = 24  ;  % surface of the scatter plot 
            in.caxis=[0 10]; % for the SNR plot 
            in.pos  = [0.13 0.11 0.7750 0.8150]; %default axis position 
            in.colorbar = true ; %add a colorbar ?
            in = lang.parse_options(in,varargin); 
            
            % energy diagram plot : 
            I0 = round(numel(h.t)/2);
            %env_matrix_t=h.t(I0:end-1);
            env_matrix_t=linspace(h.t(I0),h.t(end),size(disp_.env_matrix,1)) ;
            ax1=axes('position',in.pos) ;
            pcolor(h.disp.per,disp_.dist./env_matrix_t,s2d.norm(squeeze(disp_.env_matrix(:,:,kca))))
            shading('flat')
            xlabel('period [s]','fontweight','bold'); 
            ylabel('group vel [km/s]','fontweight','bold')
            % dispersion curve plot : 
            ax2=axes('position',in.pos) ; 
            plot(h.disp.per,disp_.gv(:,kca),'--','linewidth',in.lw,'color',in.color) ;
            scatter(h.disp.per,disp_.gv(:,kca),in.s,disp_.snr(:,kca),'fill');
            caxis(in.caxis);
            if in.colorbar 
                hc=colorbar ; 
                aa = get(hc,'Position');
                set(hc,'Position',[aa(1)+0.11 aa(2) 0.01 aa(4)]);
            end
            set(ax2,'Xtick',[]);
            set(ax2,'Ytick',[]);
            % linking the two azes : 
            linkaxes([ax1,ax2])
            ax2.Visible = 'off';
            ax2.XTick = [];
            ax2.YTick = [];
            set(ax2,'YTickLabel',{}); 
            set(ax2,'XTickLabel',{}); 
            ax2.YTick = [];
            % mking it beautiful :) 
            colormap(ax1,'gray')
            colormap(ax2,'jet')
            ylim([0.5 5]); 
            xlim([h.disp.per(1) h.disp.per(end)])
        end
        function hc = disp_pa_disp_ca(h,disp_)
        % for a a single path : plot causal/acausal disp curve
            hold on 
            h.disp_p_disp(disp_,2,'color','k'); 
            h.disp_p_disp(disp_,1,'color','r'); 
            xlabel('period [s]','fontweight','bold')
            ylabel('group vel [km/s]','fontweight','bold')
            grid minor 
            hc = colorbar;
        end
        function disp_pa_disp(h,disp_,kca,varargin)
        % for a single path : plot only causal or acausal disp curve. 
            in.color='k'; 
            in.lw   = 1 ; 
            in.s    = 24 ;    % surface of the scatter plit 
            in.toto ='' ;
            in = lang.parse_options(in,varargin);
            plot(h.disp.per,disp_.gv(:,kca),'--','linewidth',in.lw,'color',in.color) ;
            scatter(h.disp.per,disp_.gv(:,kca),in.s,disp_.snr(:,kca),'fill');
            caxis([0 10])
            ylim([1 5])
            xlim([h.disp.per(1) h.disp.per(end)])
        end
    end
    methods %reading & utility methods 
        function disp_=disp_read(h,kpath,kcmp,read_matrix) 
        % read_matrix = true : read the full diagram (if it was saved!)
            if ~exist('read_matrix') 
                read_matrix = false ; 
            end
            %read a single dispersion curve 
            % kpath : relative to h.disp 
            % jpath : relative to h. 
            jpath = h.disp.I(kpath);  
            path_= h.get_path_within_h5(jpath,kcmp) ;
            disp_     = struct()                                        ; 
            disp_.gv  = h5read(h.disp.pydb_file,['/disp_ref/gv',path_]) ; 
            disp_.snr = h5read(h.disp.pydb_file,['/disp_ref/snr',path_]); 
            disp_.lat = h.disp.lat(:,kpath);
            disp_.lon = h.disp.lon(:,kpath);
            disp_.id = h.disp.id(:,kpath)  ;
            disp_.az = h.disp.az(kpath)    ;
            disp_.baz = h.disp.baz(kpath)  ;
            disp_.dist = h.disp.dist(kpath);
            if read_matrix 
                disp_.env_matrix = h5read(h.disp.pydb_file,['/disp_ref/env_matrix',path_]);
            end
        end
        function path_=get_path_within_h5(h,I,kcmp)
            path_=['/',h.id{1,I},'/',h.id{2,I},'/',h.disp.in_ml.cmp{kcmp}];
        end
        function titre=disp_get_title(h,kpath,kcmp)
            path_name = [h.disp.id{1,kpath},'-',h.disp.id{2,kpath}];
            titre=[path_name,'  ',h.cmp{kcmp},'  ',num2str(round(h.disp.dist(kpath))),' km'];
        end
    end
    methods % monding methods ;) 
        function save_mat(h)
            dispc('reading all dispersion curves','c','b')
            h.disp.gv = h5.tree_to_matrix(h.disp.pydb_file,'/disp_ref/gv')   ;
            dispc('reading all SNR','c','b')
            [h.disp.snr ~, cmp_list]= h5.tree_to_matrix(h.disp.pydb_file,'/disp_ref/snr')   ;
            dispc('saving in .mat','c','b')
            E=struct_.add_fields_except(h,'',struct());
            E.disp.cmp=cmp_list;
            out_file=[E.disp.pydb_file(1:end-3),'.mat']    ; 
            save(out_file,'-v7.3','-struct','E')              ;
        end
        function h= disp_(pydb_file)
        % pydb_file can be a relative path to a h5 ou mat file : 
        % .. C1_01_xcorr___coher/pydb_disp__47_periods
        % .. db_disp__47_periods_with_matrix.mat
            
        % should always be the first statement 
            h = h@pycorr.live(pydb_file);
            % if this is a pydb_disp*h5 file :
            if strmatch(pydb_file(end-2:end),'.h5') 
                aa =strsplit(pydb_file,'/');
                db_file_end = aa{end}; %pydb*.h5
                db_file_end = [db_file_end(3:end-2),'mat'];
                db_file=join(aa(1:end-1),'/');
                db_file= [db_file{1},'/',db_file_end];
                % create a db_disp*mat file if it does not exist :
                if ~exist(db_file,'file')
                    keyboard
                    dispc(['  disp_ : ',db_file, 'do not exist, lets create it'],'y','b')
                    h.create_db_file(pydb_file,db_file);
                end
            else 
                db_file = pydb_file ; 
            end
            % read the db_disp*mat file 
            dispc(['  disp_ loading ',db_file],'c','n')
            E=load(db_file); 
            h.disp = E .disp ;
            
            % for some reasons the field h.disp.in_ml.cmp is not in the output pydb_disp*.h5 file !
            if ~isfield(h.disp.in_ml,'cmp')
                dispc(' disp_.m : warning : h.disp.in_ml.cmp not found -> assume we have Z-Z disp. curve','g','r');
                dispc(' if it is not the case set manually h.disp.in_ml.cmp','g','r')
                h.disp.in_ml.cmp={'ZZ'}; % so set a default value 
            end
        end
        function create_db_file(h,pydb_file,db_file);
        % load only dispersion curve metadata. Disp curve themseve will be streamed directly from HD
            [~, dset_path] = h5.tree_to_matrix(pydb_file,'/disp_ref/gv',false)   ;
            h.disp.pydb_file = pydb_file ; 
            h.disp.in_ml = h5.read_group(pydb_file,'/disp_ref/in_ml')  ;
            h.disp.md_c  = h5.read_group(pydb_file,'/disp_ref/md_c')   ;
            h.disp.in_   = h5.read_group(pydb_file,'/disp_ref/in_')    ;
            h.disp.per   = double(h.disp.md_c.periods)                  ;
            % Establish the mapping between correlation metadata and computed dispersion curve :
            % useful to display the correlation corresponding to a disp curve :
            h_id=cell_.merge(h.id(1,:),h.id(2,:),'/');
            npath = numel(dset_path);
            h.disp.I =zeros(1,npath);
            h.disp.id = dset_path;
            for ipath = 1 : numel(dset_path)
                try 
                    I = strmatch(dset_path{ipath},h_id); 
                    h.disp.I(ipath)=I(1) ;
                catch	
                    keyboard
                end
            end	
            % add extra metadata :
            h.disp.dist= h.dist(h.disp.I) ;
            h.disp.az  = h.az(h.disp.I)   ; 
            h.disp.baz = h.baz(h.disp.I)  ;
            h.disp.id  = h.id(:,h.disp.I) ; 
            h.disp.lon = h.lon(:,h.disp.I);
            h.disp.lat = h.lat(:,h.disp.I);
            % save everything w/o the dispersion curve: 
            E=struct_.add_fields_except(h,'',struct());
            save(db_file,'-v7.3','-struct','E')              ;
        end
    end
end

