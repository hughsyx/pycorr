% the measurement number ipath  (e.g h.disp.gv(:,:,ipath)) corresponds 
% to the path h.id(:,h.disp(I(ipath))). So the mapping is done by h.disp.I
% acausal, causal, sym 
classdef disp_mat 
	properties 
		disp = struct() ; 
		pydb_file=''    ;
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
		function disp_p_map(h,varargin)
		end
		function disp_p2(h,kpath,kcmp)
			ica=1;
			env_matrix = h.disp_read_matrix(kpath,1);
			I0 = round(numel(h.t)/2);
			env_matrix_t=h.t(I0:end-1);
			clf; 
			hold on
			pcolor(env_matrix_t,h.disp.per,s2d.norm(squeeze(env_matrix(:,:,2)))')
			shading('flat')
			plot(squeeze(h.disp.tsw(:,ica,kpath)),h.disp.per)
		end	
		function disp_p(h,kpath,kcmp,p1,p2,varargin)
			in.type='vel' ;% 'vel' or 'time'
			in.xlim=false; % 
			in=lang.parse_options(in,varargin);

			% convenient variable
			dist = h.disp.dist(kpath);

			% determine subplot positions : 
			dx=0.07    ;
			width = 0.25; 
			offset=0.05;
			sub(1,:) = [0.1 0.35 width 0.6]              ;
			sub(2,:) = [sub(1,1)+width+offset 0.35 width 0.6];
			sub(3,:) = [sub(2,1)+width+offset 0.35 width 0.6];

			% reading the enveloppe matrix if it exists :
			if h.disp.in_.save_matrix
				env_matrix=h.disp_read_matrix(kpath,kcmp);
				I0 = round(numel(h.t)/2);
				env_matrix_t=h.t(I0:end-1);
				env_matrix_vel = h.disp.dist(kpath)./env_matrix_t; 
			end
			for ica=[1:3]
				switch in.type 
					case{'time'} %x = temps y= periods
						xaxis = env_matrix_t ;
						yaxis = h.disp.per   ; 
						C = s2d.norm(squeeze(env_matrix(:,:,ica)))';
						xaxis2 = h.disp.dist(kpath)./h.disp.gv(:,ica,kpath); % temps arrive onde surface reinterpole
						xaxis2i= h.disp.dist(kpath)./h.disp.gvi(:,ica,kpath); % temps arrive onde original
						xlabel_='time[s]'    ;
						ylabel_='periods [s]';
						if ~in.xlim
							if dist >= 100 
								in.xlim=[0 round(dist)];
							else 
								in.xlim=[0 100];
							end
						end
				end
				ax11=axes('position',sub(ica,:)) ;
				hold on 
				pcolor(xaxis,yaxis,C);
				shading('flat');
				caxis([0 1]);
				xlim(in.xlim)
				ylim([h.disp.per(1) h.disp.per(end)]);
				if ica== 2
					xlabel(xlabel_,'fontweight','bold');
					titre=h.disp_get_title(kpath,kcmp);
					title(titre,'fontsize',12,'fontweight','bold') ;
				end
				if ica==1 
					ylabel(ylabel_,'fontweight','bold')
				end
				%..
				ax12 =axes('position',sub(ica,:));
				hold on 
				time = h.disp.dist(kpath)./h.disp.gv(:,ica,kpath);
				% per vs vel = reinterpolate signal 
				plot(xaxis2,h.disp.per,'w','linewidth',2);
				% original measurements : instant freq vs instant vel
				plot(xaxis2i,h.disp.peri(:,ica,kpath),'color',[0.5 0.5 0.5],'linewidth',2)
				% original measurements wrongly interpreted at discrete periods :
				%plot(xaxis2,h.disp.peri(:,ica,kpath),'*k','linewidth',2)
				scatter(xaxis2, h.disp.per,20,h.disp.snr(:,ica,kpath),'fill');
				xlabel(xlabel_,'fontweight','bold');
				ylabel(ylabel_,'fontweight','bold');
				colormap('jet');
				caxis([0 6])
				xlim(in.xlim)
				ylim([h.disp.per(1) h.disp.per(end)]);

				%..
				linkaxes([ax11,ax12])
				ax12.Visible = 'off';
				ax12.XTick = [];
				ax12.YTick = [];
				colormap(ax11,'jet')
				colormap(ax12,'jet')

			end
			hl=colorbar('position',[0.97 0.35 0.009 0.6]);
			ha=subplot('position',[0.13 0.1 0.775 0.15]);
			h.p_cc(h.disp.I(kpath),p1,p2,'ha',ha);
			xlimit=max(100,round(dist));
			xlim([-xlimit xlimit]);
			title('');
		end
		function titre=disp_get_title(h,kpath,kcmp)
			jpath=h.disp.I(kpath)                        ;
			path_name = [h.id{1,jpath},'-',h.id{2,jpath}];
			titre=[path_name,'  ',h.cmp{kcmp},'  ',num2str(round(h.dist(jpath))),' km'];
		end
		function env_matrix=disp_read_matrix(h,kpath,kcmp)
			jpath = h.disp.I(kpath)                                           ; 
			path_= h.get_path_within_h5(jpath,kcmp)                           ;
			env_matrix=h5read(h.disp.pydb_file,['/disp_ref/env_matrix',path_]);
		end
		function save_mat(h)
			E=struct_.add_fields_except(h,'',struct());
			out_file=[E.pydb_file(1:end-3),'.mat']    ; 
			save(out_file,'-struct','E')              ;
		end
		function h= disp_mat(mat_file)
			E=load(mat_file);
			h=struct_.add_fields_except(E,'',h)

		end
	end
end
