% C3_live : class to plot C3 using c1_trace objects:)
%   
%   C3_live_md : - deals with C3 metadata only
%                - construct a db file containing all C3 metadata 
%                - plot that used metadata only : p_map, p_station_pair,
%                - utility functions that make reading c3/c1 easier :
%
%   C1_trace  :  class to deal with a single c1_trace 
%----------------------------------------------------------------------------
% USAGE
% h=pycorr.c3_live('C3__05_c3__05-10s__coda_10dp-1000s_normenv'); 
% h.p_map 
% h.p_cc([1 2],[5 10],[10 20],1,'p_c1',false,'color','kbrg')
% clf; h.p_cc([1 2],5,10,[1 2],'p_c1',true,'color','r','c1_color','k')
%
%--------------------------------------------------------------------------
% METHODS 
%
%  p_cc : 
%    generic routing allow to plot several c3/c1 on a single or separate axis 
%      -h.read_c3 
%      -h.get_xlimit 
%........................................................................................
%  get_xlimit : 
%    determine appropriate xlimit of a plot depending on the inter-station distance
%......................................................................................
%  read_c3  :
%    read a list of c3 for different path, p1,p2, c3_type (PP/ZZ), dtype(c1/c3,...),cmp
%    trace are filtered and optionnaly normalized by c1_trace constructor. 
%    return an array of c1_trace object. 
%      -h.read_c3_single_pair_of_station(ipath,c3_type(ic3),in.dtype,in.cmp);
%........................................................................................
%  read_c3_single_pair_of_station(h,kpath,c3_type,dtype,cmp)
%    return a structure containing a raw c3 trace + metadata, ready to be converted to a c1_trace object
%...................................................................................................
% TODO : add virtual source matrix 
% normalisation uniquement au moment des plots ? 
%
% NOTE : PP should be reversed to get the right flux such as in p_cc_single_vs

%          in.dtype = 1      ; % 0='/c1' 1='/c3' 2='/c3_daily_stacked' 

classdef c3_live < pycorr.c3_live_md 
    methods 
        function p_cc_single_vs(h,kpath,p1,p2,vs,varargin) 
            %reading the vs_mat :
            c1    = h.read_c1_vs_sta(kpath,p1,p2,'vs_sta1') ;
            c1(2) = h.read_c1_vs_sta(kpath,p1,p2,'vs_sta2') ; 
            clf;
            % map :
            ha(1)=subplot(2,2,1); 
            hold on;
            in_sta=struct();
            in_sta.color=[0.9 0.9 0.9];
            in_sta.mark='^';
            in_sta.msize=14;
            
            in_map.topo=false;
            h.p_map(struct(),in_map,in_sta);
            m_plot(h.vs.lon(vs),h.vs.lat(vs),'^k','MarkerFaceColor','y','MarkerSize',14) 
            m_plot(h.lon(1,kpath),h.lat(1,kpath),'^k','MarkerFaceColor','r','MarkerSize',14);
            m_text(h.lon(1,kpath)-0.003,h.lat(1,kpath)-0.005,'1','fontsize',16,'fontweight','bold');
            m_plot(h.lon(2,kpath),h.lat(2,kpath),'^k','MarkerFaceColor','r','MarkerSize',14);
            m_text(h.lon(2,kpath)-0.003,h.lat(2,kpath)-0.005,'2','fontsize',16,'fontweight','bold');
            % c1 with window 
            ha(2)=subplot(2,2,2);
            hold on
            plot(c1(1).t,s2d.norm(c1(1).trace(:,vs))+2,'k');
            plot(c1(2).t,s2d.norm(c1(2).trace(:,vs))+0,'b');
            hl=legend([h.vs.id{vs},' - ',h.id{1,kpath}], [h.vs.id{vs},' - ',h.id{2,kpath}]);
            title('C1 between the virtual source and the 2 receivers','fontweight','bold','fontsize',16);
            set(gca,'YTick',{});
            xlabel('C1 correlation time [s]','fontweight','bold','fontSize',14);
            grid minor; 
            xlim_ = [];
            for ista = [1 : 2]
                if ista ==1 ; offset = 2; else ; offset=0; end 
                cc = c1(ista);
                patch_x1 =[cc.I_win_pos(1,vs) cc.I_win_pos(2,vs) cc.I_win_pos(2,vs) cc.I_win_pos(1,vs)];
                patch_y =[-1 -1 1 1]+offset;
                patch(cc.t(patch_x1),patch_y,'r','edgecolor','r','facealpha',0.2,'edgealpha',0.2) ;
                patch_x2 =[cc.I_win_neg(1,vs) cc.I_win_neg(2,vs) cc.I_win_neg(2,vs) cc.I_win_neg(1,vs)];
                patch(cc.t(patch_x2),patch_y,'r','edgecolor','r','facealpha',0.2,'edgealpha',0.2) ;
                xlim_ = max(abs([xlim_ ; cc.t([patch_x1 patch_x2])]));
            end
            xlim_ = min(round(5*xlim_), max(cc.t));
            xlim([-1*xlim_ xlim_]);
            aa= get(hl,'string');
            set(hl,'string',aa(1:2)) ;
            
            % Lets plot the C3 PP and NN :
            c1=h.read_c3(kpath,p1,p2,1,'dtype',0);
            ha(3)=subplot(2,2,4); 
            hold on
            % P-P :
            cc_mat = h.read_vs_matrix(kpath,p1,p2,1,3,1)   ;
            plot(h.t,flipud(s2d.norm(cc_mat.cc(:,vs)))+2,'k','linewidth',2);
            plot(c1.t,c1.trace+2,'b-','linewidth',1);
            % N-N
            cc_mat = h.read_vs_matrix(kpath,p1,p2,2,3,1)   ;
            plot(h.t,s2d.norm(cc_mat.cc(:,vs)),'k','linewidth',2);
            plot(c1.t,c1.trace,'b-');
            %
            title(['C3 using a single source dist = ', num2str(h.dist(kpath)),' km'],'fontsize',14,'fontweight','bold');
            % legend, xlim and text :
            xlimit=h.get_xlimit(kpath);
            xlim(xlimit);
            legend('C3','C1');
            xtext = xlimit(2)*80/100; 
            text(xtext,2.2,'P-P','fontSize',14,'fontweight','bold');
            text(xtext,0.2,'N-N','fontSize',14,'fontweight','bold');
            grid minor ;
            xlabel('C3 correlation time [s]','fontsize',14,'fontweight','bold')
            set(gca,'YTick',[]);
            
            %lets zoom on the SW : 
            ha(4) = subplot(2,2,3); 
            hold on 
            % P-P :
            cc_mat = h.read_vs_matrix(kpath,p1,p2,1,3,1)   ;
            plot(h.t,flipud(s2d.norm(cc_mat.cc(:,vs)))+2,'k','linewidth',2);
            plot(c1.t,c1.trace+2,'b-','linewidth',1);
            % N-N
            cc_mat = h.read_vs_matrix(kpath,p1,p2,2,3,1)   ;
            plot(h.t,s2d.norm(cc_mat.cc(:,vs)),'k','linewidth',2);
            plot(c1.t,c1.trace,'b-');
            xlim(xlimit/2);
            title('C3 with a zoom')

            
            %%
            set(ha(1),'Position',[0.05 0.58 0.4 0.33])
            set(ha(2),'Position',[0.42 0.58 0.53 0.33])
            set(ha(3),'Position',[0.42 0.11 0.53 0.33])
            set(ha(4),'Position',[0.05 0.11 0.33 0.33])

        end
        function [cc_mat]=p_cc_matrix(h,kpath,p1,p2,c3_type,varargin)
	 in.dtype = 3   ;
         in.cmp   = 1   ;
         in.norm  =true ; 
	 p_c1     =true ; 
	 in = lang.parse_options(in,varargin);
	 %reading the vs_mat :
	 cc_mat = h.read_vs_matrix(kpath,p1,p2,c3_type,in.dtype,in.cmp);
	 if p_c1 
             c1 = h.read_c3(kpath,p1,p2,c3_type,'dtype',0);
	 end
	 if in.norm 
             cc_mat.cc = s2d.norm(cc_mat.cc);
	 end
         cc_mat.cc = cc_mat.cc';
         %
         xlimit=h.get_xlimit(kpath);
	 %plot de la matrice : 
         ha = p.split_matrix_trace;
         nsta = numel(h.vs.lon);
         subplot(ha(1))
         pcolor(h.t,[1:nsta],cc_mat.cc);
         shading('flat')
         xlim(xlimit)
         title(cc_mat.titre,'interpreter','none')
         % plot de la somme : 
         subplot(ha(2))
         c3=sum(cc_mat.cc);
         plot(h.t,c3/max(abs(c3)),'b','linewidth',2)
         if p_c1 
             hold on 
             c1.p('decoration','false','color','k','norm',true);
         end
         xlim(xlimit)
     end
     function cc_mat=read_vs_matrix(h,kpath,p1,p2,c3_type,dtype,cmp)
     % normally dtype = 3 here 
         path_info = h.get_path_info(kpath,c3_type,dtype,cmp);
         cc_mat    = struct();
         cc_mat.cc = h5read(path_info.fname,path_info.path_)    ;
         cc_mat.cc = s2d.filter(cc_mat.cc,p1,p2,h.tau);
         cc_mat.cc = cc_mat.cc;
         cc_mat.path_info = path_info;
         titre = [h.id{1,kpath},'-',h.id{2,kpath}] ;
         titre = [titre,' : ',num2str(h.dist(kpath)),'km'];
         titre = [titre,' : ',path_info.prefix,' ',path_info.cmp,'-',path_info.suffix];
         cc_mat.titre = titre;
     end
     function c3=p_cc(h,kpath,p1,p2,c3_type,varargin)
     % kpath  : vector of path indice  [1,2,3,4]
     % p1     : vector of p1           [ 5,10,20]
     % p2     : vector of p2           [10,20,40]
     % c3_type: vector of c3_type : 'PP/NN',... 
         
         in.dtype = 1      ; % 0='/c1' 1='/c3' 2='/c3_daily_stacked' 
         in.cmp   = 1      ; % ZZ_ZZ, ZZ_ZN ... 
         in.norm  = true   ; % normalize each correlation ?
         in.color= 'b'     ; % color of c3 
         in.lw   = 1       ; % linewidth
         in.p_c1  = true   ; % should we plot a c1 for each c3 ? 
         in.c1_color= 'k'  ; % if yes : c1_color 
         in.c1_lw   = 2    ; % if yes : c1_linewidth 
         in.collapse=false ; % plot all trace on the same subplot ? or split them ? 
         in.ttl_fs =12     ; % legend fontzie      
         in.ttl_leg=true   ; % should we put a legend ? 
         in=lang.parse_options(in,varargin);

         % read c3 and c1 :
         c3=h.read_c3(kpath,p1,p2,c3_type,'dtype',in.dtype,'cmp',in.cmp,'norm',in.norm);
         if in.p_c1 ; 
             c1=h.read_c3(kpath,p1,p2,1,'dtype',0);
         end
         nc3 = numel(c3);
         
         %determine the color of the plot : 
         if nc3 > 1 & numel(in.color)==1 
             in.color=char(double(in.color*ones(1,nc3)));
         end
         %determine the linewidth :
         if nc3 > 1 & numel(in.lw)==1
             in.lw= in.lw*ones(1,nc3);
         end
         %determine xlimit :
         xlimit=h.get_xlimit(kpath);
         
         %sdetermine how many supbplot we have and plit the figure : 
         if in.collapse == false 
             ha = p.split_horizontal(nc3,'dy',0); 
         else
             xx = p.split_horizontal(1,'dy',0); 
             ha=ones(1,nc3);
             ha(1:end)=xx  ;
         end 
         
         %do the plot : 
         kc1 =1; 
         for icc = 1 :nc3 
             subplot(ha(icc))
             hold on
             c3(icc).p('decoration',false,'color',in.color(icc),'lw',in.lw(icc));
             if in.p_c1 
                 c3_id=[c3(icc).id{1},'_',c3(icc).id{2},'_',num2str(c3(icc).p1),'-',num2str(c3(icc).p2)];
                 c1_id=[c1(kc1).id{1},'_',c1(kc1).id{2},'_',num2str(c1(kc1).p1),'-',num2str(c1(kc1).p2)];
                 if strcmp(c3_id,c1_id) == false 
                     kc1=kc1+1;
                 end
                 c1(kc1).p('decoration',false,'color',in.c1_color,'lw',in.c1_lw)
             end
             xlim(xlimit)
         end
         xlabel('correlation time [s]','fontweight','bold','fontsize',12);
         
         %add title/legend :
         if in.collapse==true
             legend_={1,nc3};
             for icc = 1:nc3 
                 legend_(icc)={c3(icc).title}; 
             end
             hl=legend(legend_);
             legend('boxoff') 
             set(hl,'fontweight','bold','fontsize',in.ttl_fs,'interpreter','none');
         else 
             for icc=1:nc3
                 subplot(ha(icc))
                 if in.ttl_leg==false 
                     title(c3(icc).title,'fontweight','bold','fontsize',in.ttl_fs)
                 else 
                     hl=legend({c3(icc).title},'fontweight','bold','fontsize',in.ttl_fs,'interpreter','none');
                     legend('boxoff') ;
                     set(hl,'fontweight','bold','fontsize',in.ttl_fs);
                 end	
             end
         end
     end
     %-----------------------------------------
     function xlimit=get_xlimit(h,kpath)
         dist_=max(h.dist(kpath));	
         if dist_ > 100 ;
             limit= min([h.t(end) round(2*dist_)]);
             xlimit=[-limit limit];
         else
             limit =min([h.t(end) 100]);
             xlimit=[-limit limit]; 
         end	
     end	
    end
    %--------------------------------------------------------------------
    % READING METHODS : everything to read c1/c3 trace (no matrice here!)
    %------------------------------------------------------------------
    methods 
        function cc = read_c1_vs_sta(h,kpath,p1,p2,to_sta) ;
        % to_sta = vs_sta1 or vs_sta2 : indicate which c1 we read :) 
            path_info = h.get_path_info(kpath,1,1,1) ;
            aa = strsplit(path_info.path_,'/'); 
            dset_prefix = ['/c1_vs_sta/',aa{3},'/',aa{4},'/',aa{5},'/',to_sta];
            cc = struct() ;
            cc.trace     = double(h5read(path_info.fname,[dset_prefix,'/trace'])); 
            cc.trace     = s2d.filter(cc.trace,p1,p2,h.tau);
            cc.I_win_pos = double(h5read(path_info.fname,[dset_prefix,'/I_win_pos'])); 
            cc.I_win_neg = double(h5read(path_info.fname,[dset_prefix,'/I_win_neg'])); 
            cc.tau   = h.tau         ; 
            cc.t = h.t_c1a           ; 
        end
        function cc= read_c3(h,kpath,p1,p2,c3_type,varargin)
            in.dtype = 1 ;
            in.cmp  = 1  ;
            in.norm = 1  ;
            in = lang.parse_options(in,varargin); 
            % 1. get the number of path to be read : 
            npath = numel(kpath)  ; 
            nper  = numel(p1)     ;
            nc3   = numel(c3_type);
            ntot  = npath+nper+nc3;
            % 2. read them all :) 
            kcc=1;
            for ipath = 1:npath 
                for iper =1 :nper; 
                    cp1 =p1(iper); 
                    cp2 =p2(iper); 
                    for ic3 = 1:nc3 
                        trace_ = h.read_c3_single_pair_of_station(kpath(ipath),c3_type(ic3),in.dtype,in.cmp);
                        if kcc==1 
                            cc = pycorr.c1_trace(trace_,cp1,cp2,in.norm)   ;
                        else 
                            cc(kcc)=pycorr.c1_trace(trace_,cp1,cp2,in.norm);
                        end
                        kcc=kcc+1;
                    end						
                end
            end
        end
        
        %----------------------------------------------------------------------
        function cc=read_c3_single_pair_of_station(h,kpath,c3_type,dtype,cmp)
        % read the trace in the h5 file 
            path_info = h.get_path_info(kpath,c3_type,dtype,cmp);
            trace = double(h5read(path_info.fname, path_info.path_))   ;
            %get cc metadata 
            cc=[] ; 
            cc.trace = trace         ; 
            cc.id    = h.id(:,kpath) ; 
            cc.lon   = h.lon(:,kpath); 
            cc.lat   = h.lat(:,kpath); 
            cc.dist  = h.dist(kpath) ; 
            cc.az    = h.dist(kpath) ;
            cc.baz   = h.dist(kpath) ;
            cc.tau   = h.tau         ; 
            if dtype == 0 
                cc.t = h.t_c1a         ; 
            else 
                cc.t = h.t             ;
            end
            cc.path_info =  path_info; 
        end
    end
    %---------------------------------------------------
    % constructor 
    %---------------------------------------------------
    methods
        function h=c3_live(in)
            h = h@pycorr.c3_live_md(in);
        end
    end
end