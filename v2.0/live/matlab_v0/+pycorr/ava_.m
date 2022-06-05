classdef ava_ < pycorr.live
	properties 
		ava=struct()
	end
	methods 
		function ava=ava_p(h,dist1,dist2,varargin)
			in.saturation = 1; 
			in.length     = 2;  %degree
			in.med_thr    = [0 1000]; %rem all amp > med(amp)*med_thr
			in=lang.parse_options(in,varargin); 

			ava=struct(); 
			% get only path having the right distance : 
			I_dist=find(h.ava.dist >= dist1 & h.ava.dist <=dist2);
			ava.dist=h.ava.dist(I_dist);
			ava.az  =h.ava.az(I_dist)  ;
			ava.amp =h.ava.amp(I_dist) ;
			% rm meas. > some threshold : 
			I = find(ava.amp <= median(ava.amp)*in.med_thr(2) & ava.amp >= median(ava.amp)*in.med_thr(1));
			ava.amp = ava.amp(I) ; 
			ava.az  = ava.az(I)  ; 
			ava.dist= ava.dist(I);
			%sort by azimuth 
			[~, J] = sort(ava.az); 
			ava.amp=ava.amp(J)   ; 
			ava.az =ava.az(J)    ; 
			ava.dist=ava.dist(J) ;
			% get the center of the network : 
			lon = mean(h.lon(:)); 
			lat = mean(h.lat(:));
			h.p_map(struct(),struct('topo',false)) 
			p.plot_ava(ava.amp.*sqrt(ava.dist),ava.az*pi/180,in.saturation,in.length,lon,lat);
			%
		end
	end
	methods 
		function h=ava_(pydb_file);
    		h = h@pycorr.live(pydb_file)                        ;
			h.ava.ava   = h5.tree_to_matrix(pydb_file,'/ava/ava');
			h.ava.in    = h5.read_group(pydb_file,'/ava/in_')   ;
        	h.ava.in_ml = h5.read_group(pydb_file,'/ava/in_ml') ;
       		h.ava.md_c  = h5.read_group(pydb_file,'/ava/md_c')  ;
       		%flattened measurements : 
       		h.ava.amp =[h.ava.ava(1,:) h.ava.ava(2,:)];
       		h.ava.az  =[h.az h.baz];
       		h.ava.dist=[h.dist h.dist];
    	end	
	end
end	
