%% Get full dynamic model 

function [out, aux] = get_full_dyn(...
    ptpars,...
    x,...
    lhsbase,...
    lhspars,...
    prm,... 
    ref,... 
    sel,... 
    agg,... 
    gps,... 
    hivpoints,...
    sfin)

r = prm.r; p = prm.p; i = ref.i; s = ref.s; xi = ref.xi;

[r,p] = allocate_parameters(x,r,p,xi);

% Check if parameters are within bounds
tmp = [prm.bds(:,1:length(x)); x]; tmp = diff(tmp([1,3,2],:));

if min(min(tmp)) < 0 %|| p.Dx(2) > p.Dx(1) || r.default(1) > r.default(2)
    out.llk = -Inf; aux = nan;
else
    
    
    % --- Set up the necessary models -----------------------------------------
    p.housedist=prm.hcpars;
    
    p1 = p; r1=r;
    p1.IPThiv = ptpars(1);
    pt_hiv=ptpars(2);
    pt_all=ptpars(3)*(1-prm.switchhiv);
    
    iptrates=[...
        pt_all ,pt_all ,pt_hiv ];
    
    p1.IPT=iptrates;
    
    if(numel(lhspars)>6)
        r.ptdur         =  lhsbase(1);    % Duration of protection post regimen
        r.ptregdur      =  lhsbase(2);   % Regimen duration
        p.pteffi        =  lhsbase(3);    % efficacy
        p.ptdrbarrier   =  lhsbase(4);    % DR Barrier
        p.ptforg        =  lhsbase(5);    % Regimen Forgiveness
        p.ptcomp        =  lhsbase(6);    % Regimen completion (ease of adherence
        p.potency       =  lhsbase(7);    % RPotency
        
        
        r1.ptdur         =  lhspars(1);    % Duration of protection post regimen
        r1.ptregdur      =  lhspars(2);   % Regimen duration
        p1.pteffi        =  lhspars(3);    % efficacy
        p1.ptdrbarrier   =  lhspars(4);    % DR Barrier
        p1.ptforg        =  lhspars(5);    % Regimen Forgiveness
        p1.ptcomp        =  lhspars(6);    % Regimen completion (ease of adherence
        p1.potency       =  lhspars(7);    % RPotency
        
        
        
    elseif(numel(lhspars)==6)
        r1.ptdur         =  lhspars(1);    % Duration of protection post regimen
        r1.ptregdur      =  lhspars(2);   % Regimen duration
        p1.pteffi        =  lhspars(3);    % efficacy
        p1.ptforg        =  lhspars(4);    % Regimen Forgiveness
        p1.ptcomp        =  lhspars(5);    % Regimen completion (ease of adherence
        p1.potency       =  lhspars(6);    % RPotency
    end
    
    
    r.regdur1= r.ptregdur*p.forg;
    r.regdur2= r.ptregdur*(1-p.forg);
    r.pt_resist_acq =  (12./[r.regdur1 r.regdur2])*(1-p.ptdrbarrier)./p.ptdrbarrier; %Rate of acquisition of resistance during PT
    r.ptdefault     =  (12/(r.regdur1 + r.regdur2))*(1-p.ptcomp)/p.ptcomp; % estimated time lost of protection due to completion rates
    
    r.outReg        =  12./([r.regdur1 r.regdur2]);
    r.outIPT        =  12/(r.ptdur);      % 24 months protection Rangaka-Martins etc and 3HP
    
    
    r1.regdur1= r1.ptregdur*p1.forg;
    r1.regdur2= r1.ptregdur*(1-p1.forg);
    r1.pt_resist_acq =  (12./[r1.regdur1 r1.regdur2])*(1-p1.ptdrbarrier)./p1.ptdrbarrier; %Rate of acquisition of resistance during PT
    r1.ptdefault     =  (12/(r1.regdur1 + r1.regdur2))*(1-p1.ptcomp)/p1.ptcomp; % estimated time lost of protection due to completion rates
    
    r1.outReg        =  12./([r1.regdur1 r1.regdur2]);
    r1.outIPT        =  12/(r1.ptdur);      % 24 months protection Rangaka-Martins etc and 3HP
    
    
    
    % Make models 
    M1= make_model(p1, r1, i, s, gps);
    
    % 2017
    p0 = p; r0=r;
    M0 = make_model(p0, r0, i, s, gps);
    
    
    %% --- Solve the models ----------------------------------------------------
    
    % PT
    init = sfin;
    [t, soln] = ode15s(@(t,in) goveqs_scaleup(t, in, M0, M1,[2020 2020+prm.ptscale],...
        i, s,r1, p1.growth, sel, agg,hivpoints), [2019 prm.endy] , init, odeset('NonNegative',[s.nstates]));
    
    
    if (sum(any(isnan(soln)))>0)
        out.llk = -Inf; aux = nan;
    else
        
        
        
        % --- Get the objectives --------------------------------------------------
        
        
        
        tmp = interp1(t,sum(soln(:,1:i.nstates),2),t(1):t(end));
        pop = (tmp(1:end-1)+tmp(2:end))/2;
        
        
        % Prevalence
        tmp = interp1(t,sum(soln(:,s.prevalent),2),t(1):t(end));
        pre = (tmp(1:end-1)+tmp(2:end))/2;
        prev= 1e5*(pre./pop);
        
        
        % HIV Prevalence
        tmp = interp1(t,sum(soln(:,s.hivpositive),2),t(1):t(end));
        pre = (tmp(1:end-1)+tmp(2:end))/2;
        hivprev= 1e2*(pre./pop);
        
        
        %Incidence
        tmp = diff(interp1(t,soln(:,i.aux.inc),t(1):t(end)),1);
        inc = 1e5.*([tmp(:,1) , sum(tmp(:,[2 5]),2) , tmp(:,3) , tmp(:,4)]./[pop; pop;pop;pop]');
        inc_all=inc(:,1);
        inc_mdr=inc(:,2);
        inc_tbhiv=inc(:,3);
        inc_tbhivdr=inc(:,4);
        
        
        %cases
        tmp = diff(interp1(t,soln(:,i.aux.inc),t(1):t(end)),1).*p.pop1970;
        cases_all=tmp(:,1);
        cases_mdr=tmp(:,2)+tmp(:,5);
        cases_tbhiv=tmp(:,3);
        cases_tbhivdr=tmp(:,4);
        
        %deaths
        tmp = diff(interp1(t,soln(:,i.aux.mort),t(1):t(end)),1).*p.pop1970;
        deaths_all=tmp(:,2);
        deaths_tbhiv=tmp(:,3);
        deaths_mdr=tmp(:,4);
        deaths_tbhivdr=tmp(:,5);
        
        
        
        
        %Mortality
        tmpall  = diff(interp1(t,soln(:,i.aux.mort(2)),t(1):t(end)),1);
        tmphp  = diff(interp1(t,soln(:,i.aux.mort(3)),t(1):t(end)),1);
        
        mort_tbhn=1e5.*((tmpall-tmphp)./pop);
        mort_tbhp=1e5.*(tmphp./pop);
        
        
        %         tmp  = diff(interp1(t,soln(:,i.aux.notif(1)),t(1):t(end)),1);
        %         notif=1e5*(tmp./pop);
        
        
        
        %PT
        tmp = diff(interp1(t,soln(:,i.aux.pt),t(1):t(end)),1);
        tmp=tmp*p.pop1970;
        ptallpy = tmp(:,1);
        pthivpy = tmp(:,2);
        
        
        tmp = cumsum(tmp,1);
        
        ptall = tmp(:,1);
        pthiv = tmp(:,2);
        ptnohiv = tmp(:,3);
        
        
        tmp  = diff(interp1(t,soln(:,i.aux.newart(2)),t(1):t(end)),1);
        newart = tmp';
        
        tmp  = diff(interp1(t,soln(:,i.aux.newart(2)),t(1):t(end)),1).*p.pop1970;
        ltbipy = tmp';
        %
        %PT pmo
        ptpmo_all =diff(interp1 (t,soln(:,i.aux.ptpmo(1)),t(1):t(end))).*12.*p.pop1970;
        
        ptpmo_hiv =diff(interp1 (t,soln(:,i.aux.ptpmo(2)),t(1):t(end))).*12.*p.pop1970;
        
        %Dalys
        tmp =diff(interp1 (t,soln(:,i.aux.daly),t(1):t(end))).*p.pop1970;
        yld_all=tmp(:,1)*p.weights(1);
        yld_hiv=tmp(:,2)*p.weights(1);
        yll_all=tmp(:,3)*p.lex;
        yll_hiv=tmp(:,4)*p.lex;
        
        
        %PT completion
        tmp = diff(interp1(t,soln(:,i.aux.ptcomp),t(1):t(end)),1);
        tmp=tmp*p.pop1970;
        ptcompallpy = tmp(:,1);
        ptcomphivpy = tmp(:,2);
        
        %PT forgiveness (those getting forgvenes benefits after 50% pt)
        tmp = diff(interp1(t,soln(:,i.aux.ptforg),t(1):t(end)),1);
        tmp=tmp*p.pop1970;
        ptforgallpy = tmp(:,1);
        ptforghivpy = tmp(:,2);
        
        %HLM
        
        tmp = diff(interp1(t,soln(:,i.aux.pt),t(1):t(end)),1);
        tmp=tmp*p.pop1970;
        tmp = sum(tmp,1);
        hlm = tmp;
        
        
        %Notif (All)
        tmp  = diff(interp1(t,soln(:,i.aux.notif),t(1):t(end)),1);
        notif=tmp(:,1)*p.pop1970;
        
        
        out.pop1970=p.pop1970;
        out.popu=pop;
        out.prev_all=prev;
        out.inc_all=inc_all;
        
        
        out.inc_mdr=inc_mdr;
        out.inc_tbhiv=inc_tbhiv;
        out.inc_tbhivdr=inc_tbhivdr;
        out.mort_tball=mort_tbhn+mort_tbhp;
        out.mort_tbhn=mort_tbhn;
        out.mort_tbhp=mort_tbhp;
        out.notif_all=notif;
        out.ipt_all=ptall;
        out.ipt_hiv=pthiv;
        out.ipt_nohiv=ptnohiv;
        out.ipt_allpy=ptallpy;
        out.ipt_hivpy=pthivpy;
        out.iptcomp_allpy=ptcompallpy;
        out.iptcomp_hivpy=ptcomphivpy;
        out.iptforg_allpy=ptforgallpy;
        out.iptforg_hivpy=ptforghivpy;
        out.newart=newart;
        out.ltbipy=ltbipy;
        
        
        out.ptpmo_all=ptpmo_all;
        out.ptpmo_hiv=ptpmo_hiv;
        out.yld_all=yld_all;
        out.yld_hiv=yld_hiv;
        out.yll_all=yll_all;
        out.yll_hiv=yll_hiv;
        out.cases_all=cases_all;
        out.cases_mdr=cases_mdr;
        out.cases_tbhiv=cases_tbhiv;
        out.cases_tbhivdr=cases_tbhivdr;
        out.deaths_all=deaths_all;
        out.deaths_mdr=deaths_mdr;
        out.deaths_tbhiv=deaths_tbhiv;
        out.deaths_tbhivdr=deaths_tbhivdr;
        out.hivprev=hivprev;
        out.sfin=soln(end-2,1:i.nx);
        out.soln   =[soln, t];
        out.HLM =  hlm;
        
        
    end
end
