classdef live_cvg < handle 
    methods 
        function out=cvg_cc(h,kpath,varargin) 
            in = struct() ;
            in.p1=1;
            in.p2=3;
            in.kcmp =1; 
            in.ref_I=[1:numel(h.date1)]; % par defaut la ref= somme de toutes les cc
            in.cc_I =[1:numel(h.date1)]; % par defaut les cc sont prises a toutes les dates
            in.nstack =[1 10]          ; % test les cc des cc moyennees sur 1 et 10 jours 
            in.av =10                  ; % refais les mesures en tirant au sort 10 stack 
            in.time1=[0  50]           ; % interval de temps sur lequel on fait la mesure 
            in.time2=[10 100]          ; 
            in.norm = true             ; % normalize each cc prior to stacking them 
            in.norm_ref = true         ; % normalize each cc prior to define the ref 
            in.tag=''                  ; 
            in = lang.parse_options(in,varargin);
            
            % define the referece
            [cc_matrix, ref]=h.read_cc_matrix(kpath,in.p1,in.p2,in.kcmp);
            if in.norm_ref 
                ref= sum(s2d.norm(cc_matrix(:,in.ref_I)),2,'omitnan');
            else
                ref= sum(cc_matrix(:,in.ref_I),2,'omitnan');
            end
            
            % define indice of the time window 
            nwin = numel(in.time1); 
            I1 = zeros(1,nwin); 
            I2 = zeros(1,nwin);
            for iwin = 1:nwin
                I1(iwin) = find(h.t>= in.time1(iwin),1,'first');
                I2(iwin) = find(h.t<= in.time2(iwin),1,'last');
            end
            
            % loop on nstack/av/win 
            nstack = numel(in.nstack); 
            nav    = in.av;
            nI     = numel(in.cc_I);
            
            % define output structure 
            out = struct(); 
            out.ref = ref ; 
            out.cc = zeros(numel(h.t),nstack); 
            out.t  = h.t;
            out.date1= h.date1; 
            out.date2= h.date2;
            out.sta = struct();
            out.sta.lon = h.lon(:,kpath) ;
            out.sta.lat = h.lat(:,kpath); 
            out.sta.dist= h.dist(kpath); 
            out.sta.az  = h.az(kpath); 
            out.sta.baz = h.baz(kpath); 
            out.sta.id  = h.id(:,kpath);
            out.coher = zeros(nav,nstack,nwin);
            out.coher_str='nav/nstack/nwin';
            out.in = in ;
            out.name = ['+cvg/cvg_cc__',str.num2str(in.p1,3),'_',str.num2str(in.p2,3),'s'];
            out.name = [out.name,'__',num2str(out.sta.dist),'_km_',out.sta.id{1},'_',out.sta.id{2}];
            if numel(in.tag)>1;
                out.name = [out.name,'__',in.tag];
            end
            out.name = [out.name,'.mat'];
            os.mkdir('+cvg')
            for istack = 1 : nstack 
                dispc(['  working on stack #',num2str(istack),' out of ',num2str(nstack)],'c','n');
                for iav = 1 : nav 
                    %build the cc stack 
                    I_rnd = randperm(nI);
                    I_rnd = I_rnd(1:in.nstack(istack)); %iav); 
                    if in.norm
                        cc_current= sum(s2d.norm(cc_matrix(:,in.cc_I(I_rnd))),2,'omitnan');
                    else
                        cc_current= sum(cc_matrix(:,in.cc_I(I_rnd)),2,'omitnan');                       
                    end
                    if iav == 1 
                        in.cc(:,istack)=cc_current ; 
                    end
                    % compute the cc for each time window 
                    for iwin = 1 : nwin 
                        xx=corrcoef(ref(I1(iwin):I2(iwin)),cc_current(I1(iwin):I2(iwin)));
                        out.coher(iav,istack,iwin)=xx(2); 
                    end
                end
            end
            out.coher_mean = squeeze(mean(out.coher,'omitnan'));
            out.coher_std  = squeeze(std(out.coher,'omitnan'));
            save(out.name,'-struct','out');
        end
        function pws(h,kpath,p1,p2,kcmp,varargin)
            in = struct()                        ;
            in.nstack = 20                       ; 
            in.av= 10                            ; 
            in.timegate =  p2*3                  ;
            in.power =  2                        ;
            in = lang.parse_options(in,varargin) ;
            cc=h.read_cc_matrix(kpath,p1,p2,kcmp);
            [cc2 I J]= s2d.remove_zero_traces(cc);
            [nwf ntr] = size(cc2)                ;
            % dbg 
            max_tr = 300 ; 
            ntr = min(max_tr,ntr);
            cc2 = cc2(:,[1:min(max_tr,ntr)]);
            %
            ref = sum(cc2,2)                     ;
            %prestack
            D = zeros(nwf,in.nstack); 
            for istack = 1 : in.nstack 
                for iav = 1 : in.av 
                    I_tr = randperm(ntr,istack)                 ;
                    D(:,istack) = D(:,istack)+sum(cc2(:,I_tr),2);
                end
            end
            % stacking the inst. phase of the signal : 	
            Da=angle(hilbert(D))          ;
            phase_stack.real= zeros(1,nwf); 
            phase_stack.imag= zeros(1,nwf); 
            for istack= 1: in.nstack
                phase_stack.real = phase_stack.real + cos(Da(:,istack))';
                phase_stack.imag = phase_stack.imag + sin(Da(:,istack))';
            end
            % determine the coherency :
            coh = 1. / in.nstack * abs(complex(phase_stack.real,phase_stack.imag));
            timegate_npts = round(in.timegate * 1/h.tau);
            boxcar=rectwin(timegate_npts)/timegate_npts;
            coh2=conv(coh,boxcar,'same');
            coh2 = coh2.^in.power ;
            % getting new refs : 
            ref_pws = zeros(nwf,1);
            for istack = 1 : in.nstack
                ref_pws = ref_pws + (D(:,istack).* coh2');
            end
            whos cc2
            clf; hold on ; plot(h.t,s2d.norm(ref_pws),'b','linewidth',2) ; plot(h.t,s2d.norm(ref),'k')
            xlim([-1000 1000])
        end
        function [snr success]=cvg_(h,kpath,p1,p2,kcmp,varargin)
            in.stackv=[]; % vecteur contenant le nb jour a stacker
            in.av = 10  ; 
            in = lang.parse_options(in,varargin) ;
            %getting default argument :
            if isempty(in.stackv)
                nday = numel(h.date1);
                in.stackv = [2:nday];
            end
            
            % read data : 
            cc=h.read_cc_matrix(kpath,p1,p2,kcmp);
            [cc2 I J]= s2d.remove_zero_traces(cc);
            ref = sum(cc2,2)        ;
            [nwf ntr] = size(cc2)   ;
            I=find(in.stackv<= ntr) ; 
            in.stackv = in.stackv(I);
            
            % init output : 
            nstack = numel(in.stackv);
            snr.snr =zeros(nstack,2);  %acausal - causal
            snr.in.kpath = kpath ;
            snr.in.kcmp  = kcmp  ;
            snr.in.p1    = p1    ; 
            snr.in.p2    = p2    ; 
            snr.in.stackv= in.stackv  ;
            snr.in.av    = in.av      ;
            snr.lon   = h.lon(:,kpath);
            snr.lat   = h.lat(:,kpath); 
            snr.id    = h.id(:,kpath) ; 
            snr.dist  = h.dist(kpath) ; 
            snr.az    = h.az(kpath)   ;
            snr.baz   = h.baz(kpath)  ; 
            snr.ref   = ref ; 
            snr.t     = h.t ; 
            titre=[snr.id{1},'-',snr.id{2},'   ',str.num2str(snr.dist,4),'km']; 
            titre=[titre,'   ',str.num2str(p1,3),'-',str.num2str(p2,3),'s'];
            titre=[titre,'   kcmp=',num2str(kcmp),' ',h.cmp{kcmp} ] ;
            snr.titre=titre ;
            
            % compute snr vs number of day : 
            kstack=0;
            for istack = in.stackv 
                kstack = kstack + 1    ;
                D    = zeros(nwf,in.av);
                for iav = 1 : in.av 
                    I_tr = randperm(ntr,istack)  ;
                    D(:,iav) = sum(cc2(:,I_tr),2);
                end
                snr_D=cc2d.sw_snr(D,ref,h.t,h.dist(kpath),'p1',p1,'p2',p2);
                snr.snr(kstack,:) = median(log(snr_D.amp./snr_D.noiz),2)  ;
            end
            if exist('snr_D')
                snr.Isurf1 = snr_D.Isurf1(:,1) ; 
                snr.Isurf2 = snr_D.Isurf2(:,1) ; 
                snr.Inoise1 = snr_D.Inoise1(:,1);
                snr.Inoise2 = snr_D.Inoise2(:,1);
                snr.I1 = snr_D.I1(:,1) ; 
                snr.I2 = snr_D.I2(:,1) ;
                success = true;
            else
                success =false ;
            end
        end
    end
end