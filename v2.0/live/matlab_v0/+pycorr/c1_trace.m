classdef c1_trace < handle 
	properties
		trace = [] ;
		id    = [] ;
		lon   = [] ; 
		lat   = [] ; 
		dist  = [] ; 
		az    = [] ; 
		baz   = [] ;
    tau   = [] ;
		t     = [] ;
		title = '' ; 
		p1    = [] ; 
		p2    = [] ; 
		cmp   = '' ;
		path_info=struct();
	end
	methods
		function p_fft(h,varargin)
		end
		function p(h,varargin)
			in.decoration = true ;
			in.color = 'k'       ;
			in.lw    = 1         ;
			in.norm  = true      ;
			in = lang.parse_options(in,varargin); 

			if in.norm 
				trace_=h.trace./max(abs(h.trace)); 
			else 
				trace_=h.trace; 
			end
			plot(h.t,trace_,'color',in.color,'linewidth',in.lw)
			if in.decoration 
				p.label('correlation time [s]','',h.title)
			end
		end
	end
	methods 
		function filter(h,p1,p2); 
			h.trace = s2d.filter(h.trace,p1,p2,h.tau);
		end
		function mk_title(h); 
			titre= [h.id{1},'-',h.id{2}];
			titre=[titre,' : ',num2str(round(h.dist)),'km'];
			titre=[titre,' : ',h.path_info.prefix(2:end)];
			titre=[titre,'-',h.path_info.cmp];
			if ~isempty(h.path_info.suffix)
				titre= [titre,'-',h.path_info.suffix];
			end
			titre=[titre,'  ',num2str(h.p1),'-',num2str(h.p2),'s'];
			h.title=titre;
		end
			function h=c1_trace(cc,p1,p2,norm)
			h=struct_.copy_field_to(cc,h);
			h.p1  = p1;
			h.p2  = p2; 
			h.filter(p1,p2); 
			if norm 
				h.trace = s2d.norm(h.trace);
			end
			h.mk_title
		end
	end
end