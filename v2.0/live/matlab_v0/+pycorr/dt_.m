% dt_measure : measure median/mean dt for a single station for a given range of dist, and corr coeff.
%              return a dt structure containing all the information 
%              to count the number of measurements : assume that a measurement was not done if corr coeff==0.
classdef dt_ < pycorr.live 
	properties 
		dt = struct()
	end
	methods
		function dt = p_dt_single_station(h,ksta,varargin)
			in.p_median=true; 
      in.median_color='k'; 
      in.median_lw=2;
      in.type='scatter'; %or 'proba' 
      in.scatter_size=4; 
      in.ylim=[];
      in.xlim=[];
      in.coeff=[-1 1] ;
      in.dist =[0.1 inf]; 
      in.p_nmeas = true ;
      in=lang.parse_options(in,varargin); 

      nsubplot=1;
      if in.p_nmeas == true 
      	nsubplot =2 ;
      end
      if nsubplot==1
	      ha=p.split_horizontal(nsubplot,'dy',0) ;
	    else 
	    	ha = p.split_matrix_trace; 
	    end
      dt   =h.dt_measure(ksta,'dist',in.dist,'coeff',in.coeff);
      nsta =size(dt.matrix,2);
      xlimit=[dt.date_(1) dt.date_(end)];

      subplot(ha(1))
      hold on 
      switch in.type 
	      case{'scatter'}
	      	for ista=1:nsta 
	      		I=find(dt.coeff_matrix(:,ista) >0);
  	    		scatter(dt.date_(I),dt.matrix(I,ista),in.scatter_size,dt.coeff_matrix(I,ista),'filled');
    	    end
      	  colorbar
       	 caxis([0 1])
       	 if ~isempty(in.ylim); 
       	 	ylim(in.ylim)
       	 end
       	 if isempty(in.xlim)
       	 	xlim(xlimit)
       	 else 
       	 	xlim(in.xlim)
       	 end
      end
      grid minor
      datetick('x','mmm','keeplimits','keepticks');
      p.label('calendar time','dt [s]',dt.title) 

      if in.p_median 
        plot(dt.date_,dt.median,'color',in.median_color,'linewidth',in.median_lw);
      end
      %.....................
      if in.p_nmeas == true 
      	xlabel('');
      	subplot(ha(2))
      	plot(dt.date_,dt.nmeas,'color',in.median_color,'linewidth',2);
      	if isempty(in.xlim)
       	 	xlim(xlimit)
       	else 
       	 	xlim(in.xlim)
       	end
      	hc=colorbar            ; %just to sync both plot;
      	set(hc,'visible','off'); 
       	grid minor 
      	datetick('x','mmm','keeplimits','keepticks')
      	p.label('calendar time','nmeas','');
      end
		end
	end
	methods 
		function dt = dt_measure(h,ksta,varargin)
	    in.coeff=[-inf inf]; %select path having this range of coeff
      in.dist =[0.1 inf]   ;
      in=lang.parse_options(in,varargin);
      dt.in = in ;

      [I1 I2 sta_info]=h.get_indice_path_with_station(ksta);
      I=[I1 I2]                                   ; 
      dt_ksta_net = h.dt.dt(:,I1)*-1              ;
      dt_net_ksta = h.dt.dt(:,I2)                 ;
      ndate = size(h.dt.coeff,1)                  ; 
      nsta  = numel(h.get_station_list.lon)       ;
      dt.matrix_no_selection = [dt_ksta_net dt_net_ksta];
      dt.coeff_matrix_no_selection = h.dt.coeff(:,I);
      dt.matrix = zeros(ndate,nsta)               ;
      dt.coeff_matrix = zeros(ndate,nsta)         ;
      dt.date_=h.dt.date                          ; 
      dt.sta = sta_info                           ;
      % set additional title string : 
      dt.title ={} ;
      dt.title(1) = {['clock drift for ',dt.sta.name]};
      dt.title(2) = {['coeff = [',num2str(in.coeff(1)),'-',num2str(in.coeff(2)),']']};
      dt.title(2)  = {[dt.title{2},'    dist = [',num2str(in.dist(1)),'-',num2str(in.dist(2)),']']};
      % add distance between ksta and the rest of the network : 
      dt.dist = h.dist(I) ;
      %scan each day of measurements to select appropriate paths :
      nday=size(dt.matrix,1)         ;
      dt.median       = zeros(1,nday);
      dt.mean         = zeros(1,nday);
      dt.coeff_median = zeros(1,nday);
      dt.coeff_mean   = zeros(1,nday);
      dt.nmeas        = zeros(1,nday);
      I_dist  = find(dt.dist >= in.dist(1) & dt.dist <= in.dist(2));
      for iday = [1 :nday]
      	coeff_day = dt.coeff_matrix_no_selection(iday,:)           ;
      	I_coeff = find(coeff_day >= in.coeff(1) & coeff_day <= in.coeff(2));
      	I_sel = intersect(I_dist, I_coeff)                         ;
      	dt.coeff_matrix(iday,I_sel) = coeff_day(I_sel)             ;
      	dt.matrix(iday,I_sel)       = dt.matrix_no_selection(iday,I_sel);
      	dt.median(iday) = median(dt.matrix(iday,I_sel))            ;
        dt.mean(iday)   = mean(dt.matrix(iday,I_sel))              ;
        dt.coeff_mean(iday)   = mean(dt.coeff_matrix(iday,I_sel))  ;
        dt.coeff_median(iday) = median(dt.coeff_matrix(iday,I_sel));
        dt.nmeas(iday) = numel(find(coeff_day(I_sel) ~= 0))        ;
      end
    end
		%function dt_median=dt_get_median_value_date_interval(h,ksta,varargin)
    %  in.coeff=[0 1]; 
    %  in.dist=[0 inf];
    %  in.date1=datenum(2012,1,1);
    %  in.date2=datenum(2013,1,1);
    %  in=lang.parse_options(in,varargin);
    %  dt=h.dt_measure(ksta,'coeff',in.coeff','dist',in.dist);
    %  
    %  Idate=find(dt.date_ >=in.date1 & dt.date_ <= in.date2);
    %  Inum =find(isnan(dt.median)==0);
    %  dt_median=median(dt.median(intersect(Idate,Inum)));
    %end

	end
	methods 
		function h=dt_(pydb_file)
    		h = h@pycorr.live(pydb_file)                         ;
    		h.dt.dt    = h5.tree_to_matrix(pydb_file,'/dt/dt')   ;
        h.dt.coeff = h5.tree_to_matrix(pydb_file,'/dt/coeff');
        h.dt.in    = h5.read_group(pydb_file,'/dt/in_')      ;
        h.dt.in_ml = h5.read_group(pydb_file,'/dt/in_ml')    ;
        h.dt.md_c  = h5.read_group(pydb_file,'/dt/md_c')     ;
        h.dt.date1 = h.dt.md_c.date1+366 ;                                
        h.dt.date2 = h.dt.md_c.date2+366 ;
        h.dt.date  = h.dt.md_c.date +366 ;
        %tmp fix ; 
        h.dt.dt =squeeze(h.dt.dt);
        h.dt.coeff= squeeze(h.dt.coeff);
    	end
	end
end	