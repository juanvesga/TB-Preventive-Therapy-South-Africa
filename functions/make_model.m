function M = make_model(p, r, i, s, gps)

% --- Get the linear rates ------------------------------------------------
m = zeros(i.nstates);
mh = zeros(i.nstates);
mc = zeros(i.nstates);
mfo=zeros(i.nstates);
 


%---- TB cascade


for ig = 1:length(gps.hiv)
    hiv = gps.hiv{ig};
    
    gu = @(st) i.(st).(hiv);
    U     = gu('U');
    Pu1   = gu('Pu1');
    Pu2   = gu('Pu2');
    Qu    = gu('Qu');
    Ru    = gu('Ru');
    Qc    = gu('Qc');

    
    
    
    
    % IPT
    source = U; destin = Pu1; rate = (1-p.LTBItest)*p.IPT(ig)*p.housedist(1);
    m(destin, source) = m(destin, source) + rate;

    % Regimen completion (no missed doses)
    source = Pu1; destin = Pu2; rate = r.outReg(1);
    m(destin, source) = m(destin, source) + rate;

    
    
    % Regimen completion (no missed doses)
    source = Pu2; destin = Qu; rate = r.outReg(2);
    m(destin, source) = m(destin, source) + rate;
    mc(destin, source) = mc(destin, source) + rate;
    
    
    % End of PT protection
    source = Qu; destin = Ru; rate = r.outIPT;
    m(destin, source) = m(destin, source) + rate;
    
    
    %PT default (during PT administration) no forgiveness
    source = Pu1; destin = Ru; rate = r.ptdefault;
    m(destin, source) = m(destin, source) + rate;
    
    source = Pu2; destin = Ru; rate = r.ptdefault*(1-p.ptforg);
    m(destin, source) = m(destin, source) + rate;
    
    %PT default (during PT administration) with forgiveness
    source = Pu2; destin = Qu; rate = r.ptdefault*p.ptforg;
    m(destin, source) = m(destin, source) + rate;
    mfo(destin, source) = mfo(destin, source) + rate;
    
    
    
    for istr = 1:length(gps.strain)
        strain = gps.strain{istr};
        
        gi = @(st) i.(st).(strain).(hiv);
        Lf    = gi('Lf');
        Ls    = gi('Ls');
        Pf1    = gi('Pf1');
        Ps1    = gi('Ps1');
        Pf2    = gi('Pf2');
        Ps2    = gi('Ps2');
        Qf    = gi('Qf');
        Qs    = gi('Qs');
        Rf    = gi('Rf');
        Rs    = gi('Rs');
        Ia    = gi('Ia');
        Is    = gi('Is');
        Dx    = gi('Dx');
        Tx    = gi('Tx');
        Tx2   = gi('Tx2');
        E     = gi('E');
        Rlo   = gi('Rlo');
        Rhi   = gi('Rhi');
        R     = gi('R');
        
        % Group Selectors
        ishiv = (ig > 1);
        ishivcas = (ig > 2);
        ishivart = (ig > 2);
        ismdr   = strcmp(strain, 'mdr');
        
        
        % Recruitment into PT
        source = Lf; destin = Pf1; rate = p.IPT(ig)*p.housedist(2);
        m(destin, source) = m(destin, source) + rate;
        
        source = Ls; destin = Ps1; rate = p.IPT(ig)*p.housedist(3);
        m(destin, source) = m(destin, source) + rate;
        
        % Regimen completion
        source = Pf1; destin = Pf2; rate = r.outReg(1);
        m(destin, source) = m(destin, source) + rate;
        
        source = Ps1; destin = Ps2; rate = r.outReg(1);
        m(destin, source) = m(destin, source) + rate;

        source = Pf2; destin = Qf; rate = r.outReg(2);
        m(destin, source) = m(destin, source) + rate;
        mc(destin, source) = mc(destin, source) + rate;
        
        
        source = Ps2; destin = Qs; rate = r.outReg(2);
        m(destin, source) = m(destin, source) + rate;
        mc(destin, source) = mc(destin, source) + rate;
        
        % Regimen cure by regimens potency
        if (p.potency>0)
        rate= min(52,((12/(r.regdur1 + r.regdur2))*p.potency/(1-p.potency))*(ismdr~=1));      
        else
        rate=0;
        end
        
        source = Pf1; destin = Qc; 
        m(destin, source) = m(destin, source) + rate*r.pot_split;
        mc(destin, source) = mc(destin, source) + rate*r.pot_split; 
         
        source = Pf2; destin = Qc; 
        m(destin, source) = m(destin, source) + rate;
        mc(destin, source) = mc(destin, source) + rate;
        
        source = Ps1; destin = Qc; 
        m(destin, source) = m(destin, source) + rate*r.pot_split;
        mc(destin, source) = mc(destin, source) + rate*r.pot_split;
        
        source = Ps2; destin = Qc; 
        m(destin, source) = m(destin, source) + rate;
        mc(destin, source) = mc(destin, source) + rate;
        
        
        % End of PT protection
        source = Qf; destin = Rf; rate =r. outIPT;
        m(destin, source) = m(destin, source) + rate;
        
        source = Qs; destin = Rs; rate = r.outIPT;
        m(destin, source) = m(destin, source) + rate;
        
        %DR acquisition while on PT
%         source = Pf; destin = i.('Pf').('mdr').(hiv);
%         rate = r.pt_resist_acq*(ismdr~=1);
%         m(destin, source) = m(destin, source) + rate;
%         
%         source = Ps; destin = i.('Ps').('mdr').(hiv);
%         rate = r.pt_resist_acq*(ismdr~=1);
%         m(destin, source) = m(destin, source) + rate;
        
        %PT default from first half of regimen
        source = Pf1; destin = Rf; 
        rate = r.ptdefault*p.ptdrbarrier;
        m(destin, source) = m(destin, source) + rate;
        
        source = Pf1; destin = i.('Rf').('mdr').(hiv); 
        rate = r.ptdefault*(1-p.ptdrbarrier);
        m(destin, source) = m(destin, source) + rate;

        source = Ps1; destin = Rf; 
        rate = r.ptdefault*p.ptdrbarrier;
        m(destin, source) = m(destin, source) + rate;
        
        source = Ps1; destin = i.('Rf').('mdr').(hiv); 
        rate = r.ptdefault*(1-p.ptdrbarrier);
        m(destin, source) = m(destin, source) + rate;


        %PT default from 2nd half (during PT administration) no forgiveness
        source = Pf2; destin = Rf; 
        rate = r.ptdefault*(1-p.ptforg)*p.ptdrbarrier;
        m(destin, source) = m(destin, source) + rate;
        
        source = Pf2; destin = i.('Rf').('mdr').(hiv); 
        rate = r.ptdefault*(1-p.ptforg)*(1-p.ptdrbarrier);
        m(destin, source) = m(destin, source) + rate;
        
        
        source = Ps2; destin = Rs; 
        rate = r.ptdefault*(1-p.ptforg)*p.ptdrbarrier;
        m(destin, source) = m(destin, source) + rate;
        
        source = Ps2; destin = i.('Rs').('mdr').(hiv); 
        rate = r.ptdefault*(1-p.ptforg)*(1-p.ptdrbarrier);
        m(destin, source) = m(destin, source) + rate;
 
        %PT default from 2nd half (during PT administration) with forgiveness No cure
        source = Pf2; destin = Qf; 
        rate = r.ptdefault*(p.ptforg)*p.ptdrbarrier;
        m(destin, source) = m(destin, source) + rate;
        mfo(destin, source) = mfo(destin, source) + rate;
                
        source = Pf2; destin = i.('Qf').('mdr').(hiv); 
        rate = r.ptdefault*(p.ptforg)*(1-p.ptdrbarrier);
        m(destin, source) = m(destin, source) + rate;
        mfo(destin, source) = mfo(destin, source) + rate;
        
        source = Ps2; destin = Qs; 
        rate = r.ptdefault*(p.ptforg)*p.ptdrbarrier;
        m(destin, source) = m(destin, source) + rate;
        mfo(destin, source) = mfo(destin, source) + rate;
        
        source = Ps2; destin = i.('Qs').('mdr').(hiv); 
        rate = r.ptdefault*(p.ptforg)*(1-p.ptdrbarrier);
        m(destin, source) = m(destin, source) + rate;
        mfo(destin, source) = mfo(destin, source) + rate;
        
        % Fast to slow track
        source = Lf; destin = Ls; rate = r.slow;%*(1-ishiv);
        m(destin, source) = m(destin, source) + rate;
        
        source = Pf1; destin = Ps1; rate = r.slow;%*(1-ishiv);
        m(destin, source) = m(destin, source) + rate;
        
        source = Pf2; destin = Ps2; rate = r.slow;%*(1-ishiv);
        m(destin, source) = m(destin, source) + rate;
        
        source = Qf; destin = Qs; rate = r.slow;%*(1-ishiv);
        m(destin, source) = m(destin, source) + rate;
        
        source = Rf; destin = Rs; rate = r.slow;%*(1-ishiv);
        m(destin, source) = m(destin, source) + rate;
        
        
        % --- Reactivation No PT effect
        RRh = max(1, r.progRRhiv*ishiv)*max(1-ishivart, r.ARTred); % Relative reduction for HIV
        
        source = Lf; destin = Ia; rate = r.fast_react*RRh;
        m(destin, source) = m(destin, source) + rate;
        
        
        source = Ls; destin = Ia; rate = r.slow_react*RRh;
        m(destin, source) = m(destin, source) + rate;
        
        source = Rf; destin = Ia; rate = r.fast_react*RRh;
        m(destin, source) = m(destin, source) + rate;
        
        source = Rs; destin = Ia; rate = r.slow_react*RRh;
        m(destin, source) = m(destin, source) + rate;
        
        
        % --- Reactivation with PT effect
        suppress=ishiv*min(p.pteffi/r.corr_hiv,1) + p.pteffi*(1-ishiv);
        source = Pf1; destin = Ia;
        rate = r.fast_react*RRh*ismdr*(1-suppress*p.effondr) + r.fast_react*RRh*(1-ismdr)*(1-suppress);
        m(destin, source) = m(destin, source) + rate;
        mc(destin, source) = mc(destin, source) + rate;
        
        source = Ps1; destin = Ia;
        rate = r.slow_react*RRh*ismdr*(1-suppress*p.effondr) + r.slow_react*RRh*(1-ismdr)*(1-suppress);
        m(destin, source) = m(destin, source) + rate;
        mc(destin, source) = mc(destin, source) + rate;

        source = Pf2; destin = Ia;
        rate = r.fast_react*RRh*ismdr*(1-suppress*p.effondr) + r.fast_react*RRh*(1-ismdr)*(1-suppress);
        m(destin, source) = m(destin, source) + rate;
        mc(destin, source) = mc(destin, source) + rate;
        
        source = Ps2; destin = Ia;
        rate = r.slow_react*RRh*ismdr*(1-suppress*p.effondr) + r.slow_react*RRh*(1-ismdr)*(1-suppress);
        m(destin, source) = m(destin, source) + rate;
        mc(destin, source) = mc(destin, source) + rate;

        
        
        source = Qf; destin = Ia;
        rate = r.fast_react*RRh*ismdr*(1-suppress*p.effondr) + r.fast_react*RRh*(1-ismdr)*(1-suppress);
        m(destin, source) = m(destin, source) + rate;
        
        source = Qs; destin = Ia;
        rate = r.slow_react*RRh*ismdr*(1-suppress*p.effondr) + r.slow_react*RRh*(1-ismdr)*(1-suppress);
        m(destin, source) = m(destin, source) + rate;
        
        
        % Breakdown Ia to Is
        source = Ia; destin = Is; rate = r.symp_del;
        m(destin,source) = m(destin,source) + rate;
        
        
        % --- Primary careseeking, including access to public sector care
        source = Is; destin = Dx;
        rate = r.careseek*max(1, r.csRRhiv*ishivcas);
        csr=rate;
        m(destin,source) = m(destin,source) + rate;
        
        % -- Get transitions by providers
        
        DxAttempt = p.Dx(ig);
        
        pFLinit= DxAttempt*p.Tx_init*( ((1-ismdr)*p.xpert*p.xpert_sens) + ((1-p.xpert)*p.smear_sens) );
        
        pSLinit= DxAttempt*p.Tx_init*ismdr*p.xpert*p.xpert_sens;
        
        p_ltfu  = 1-(DxAttempt*p.Tx_init* ((p.xpert*p.xpert_sens) + (1-p.xpert)*p.smear_sens));
        
        source  = Dx;
        destins =       [Tx,       Tx2,   E];
        rates   = r.Dx* [pFLinit,  pSLinit, p_ltfu];
        m(destins, source) = m(destins, source) + rates';
        
        
        % --- FL Treatment
        
        pFLcure = p.cure;
        rMDRacq = r.MDR_acqu;
        pSLtran = p.SL_trans;
        
        if ismdr~=1
            
            source  = Tx;
            destins = [Rlo            E                  Rhi,                     i.Tx.mdr.(hiv)];
            rates   = [r.Tx*pFLcure,  r.Tx*(1-pFLcure),  r.default*(1-rMDRacq),   r.default*rMDRacq];
            m(destins, source) = m(destins, source) + rates';
            
        elseif  ismdr
            
            source  = Tx;
            destins = [Tx2             E];
            rates   = [r.Tx*pSLtran,   r.Tx*(1-pSLtran) + r.default];
            m(destins, source) = m(destins, source) + rates';
            
        end
        
        % --- SL Treatment
        source  = Tx2;
        destins = [Rlo                 E                     ];
        rates   = [r.Tx2*p.cure2,     r.Tx2*(1-p.cure2)+  r.default2];
        m(destins, source) = m(destins, source) + rates';
        
        
        % --- Secondary careseeking
        source = E; destin = Dx;
        rate =ishivcas*csr  + (1-ishivcas)*csr*r.cs2;
        m(destin, source) = m(destin, source) + rate;
        
        % --- Relapse
        sources = [Rlo, Rhi, R];
        destin  = Ia;
        rates   = r.relapse;
        m(destin, sources) = m(destin, sources) + rates;
        
        sources = [Rlo, Rhi];
        destin  = R;
        rates   = 0.5;
        m(destin, sources) = m(destin, sources) + rates;
        
        
        % --- Self cure
        rate = r.selfcure;
        source = [Ia, Is,E,Dx];%intersect(s.infectious,s.(age));
        destin = Rlo;
        m(destin, source) = m(destin, source) + rate;
        
        
        %             rate = r.selfcure;
        %             for source = intersect(s.infectious,s.(age))
        %                 destin = Rlo;
        %                 m(destin, source) = m(destin, source) + rate;
        %             end
    end
end



M.lin   = sparse(m - diag(sum(m,1)));
M.ptcomp= sparse(mc - diag(sum(mc,1)));
M.ptforg= sparse(mfo - diag(sum(mfo,1)));


% --- Get HIV cascade linear transitions  ---------------------------------------------
%- to ART
sources = intersect(s.hu,s.U); destins = intersect(s.ht,s.U);
rates   = p.hivtest_notb*p.art_notb * p.ARToff*r.ARTrec*r.OptimART*(1-p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;

sources = intersect(s.hu,s.U); destins = intersect(s.ht,s.Pu1);
rates   = p.hivtest_notb*p.art_notb * p.ARToff*r.ARTrec*r.OptimART*(p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;

sources = intersect(s.hu,s.Lf); destins = intersect(s.ht,s.Lf);
rates   = p.hivtest_notb*p.art_notb * p.ARToff*r.ARTrec*r.OptimART*(1-p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;

sources = intersect(s.hu,s.Lf); destins = intersect(s.ht,s.Pf1);
rates   = p.hivtest_notb*p.art_notb * p.ARToff*r.ARTrec*r.OptimART*(p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;


sources = intersect(s.hu,s.Ls); destins = intersect(s.ht,s.Ls);
rates   = p.hivtest_notb*p.art_notb * p.ARToff*r.ARTrec*r.OptimART*(1-p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;

sources = intersect(s.hu,s.Ls); destins = intersect(s.ht,s.Ps1);
rates   = p.hivtest_notb*p.art_notb * p.ARToff*r.ARTrec*r.OptimART*(p.IPThiv);
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;


sources = intersect(s.hu,s.notbcare); destins = intersect(s.ht,s.notbcare);
rates   = p.hivtest_notb*p.art_notb * p.ARToff*r.ARTrec*r.OptimART;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;

sources = intersect(s.hu,s.tbcare); destins = intersect(s.ht,s.tbcare);
rates   = p.hivtest_tb*p.art_tb * p.ARToff*r.ARTrec*r.OptimART;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;

%ART dropout
sources = s.ht; destins = s.hu;
rates   = r.art_dropout;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
mh(inds) = mh(inds) + rates;

M.linh = sparse(mh - diag(sum(mh,1)));

% --- Get nonlinear rates of TB transmission
for istr = 1:length(gps.strain)
    strain = gps.strain{istr};
    
    
    m = zeros(i.nstates);
    
    for ig = 1:length(gps.hiv)       
        hiv = gps.hiv{ig};
        
        getso = @(st) [i.(st).('ds').(hiv) , i.(st).('mdr').(hiv) ];
        getdes = @(st) i.(st).(strain).(hiv);
        
        U = i.U.(hiv);
        Pu1 = i.Pu1.(hiv);
        Qu = i.Qu.(hiv);
        Ru = i.Ru.(hiv);
        
        Qc = i.Qc.(hiv);
        
        Lf  = getdes('Lf');
        Pf1  = getdes('Pf1');
        Pf2  = getdes('Pf2');
        Qf  = getdes('Qf');
        Rf  = getdes('Rf');
        
        
        Lfso  = getso('Lf');
        Pf1so = getso('Pf1');
        Pf2so = getso('Pf2');
        Qfso = getso('Qf');
        Rfso  = getso('Rf');
        Ls  = getso('Ls');
        Ps1 = getso('Ps1');
        Ps2 = getso('Ps2');
        Qs = getso('Qs');
        Rs  = getso('Rs');
        Rlo = getso('Rlo');
        Rhi = getso('Rhi');
        R   = getso('R');
        
        m(Lf, [U Lfso Ls Rlo Rhi R]) = 1;
        m(Pf1, [Pu1 Pf1so Ps1]) = 1;
        m(Pf2, [Pu2 Pf2so Ps2]) = 1;
        m(Qf, [Qu Qfso Qc Qs]) = 1;
        m(Rf, [Ru Rfso Rs]) = 1;
        
    end
    %HIV Infection susceptibilty and immunity
    m(:,[ s.Lf s.Ls s.ipt_all s.Qc s.Qf s.Qs s.Rc s.Rf s.Rs...
        s.Rlo s.Rhi s.R]) = m(:,[ s.Lf s.Ls s.ipt_all...
        s.Qc s.Qf s.Qs s.Rc s.Rf s.Rs s.Rlo s.Rhi s.R])*p.imm;
    m(:,s.hu) = m(:,s.hu)*r.reinfhiv;
    
    m(intersect(s.ds,s.hn),s.ipt_all) = m(intersect(s.ds,s.hn),s.ipt_all)*(1-p.pteffi*p.reinfprotection);
    m(intersect(s.ds,s.hivpositive),s.ipt_all) = m(intersect(s.ds,s.hivpositive),s.ipt_all)*(1-min(p.pteffi/r.corr_hiv,1)*p.reinfprotection);
    m(intersect(s.mdr,s.hn),s.ipt_all) = m(intersect(s.mdr,s.hn),s.ipt_all)*(1-p.pteffi*p.reinfprotection*p.effondr);
    m(intersect(s.mdr,s.hivpositive),s.ipt_all) = m(intersect(s.mdr,s.hivpositive),s.ipt_all)*(1-min(p.pteffi/r.corr_hiv,1)*p.effondr*p.reinfprotection);
    
    
    M.nlin.(strain) = sparse(m - diag(sum(m,1)));
end




% Non linear HIV acquisition
m = zeros(i.nstates);
sources = s.hn;
destins = s.hu;
rates   = p.hivcq_ad;
inds    = sub2ind([i.nstates, i.nstates],destins,sources);
m(inds) = m(inds) + rates;

M.nlinh = sparse(m - diag(sum(m,1)));


% --- Getting force-of-infection for children and adults settings
m = zeros(2,i.nstates);%1. Chi-Ds 2.Chi-MDR 3.Ad-DS 4.Adu-MDR
m(1,intersect(s.infectious,s.ds))  = r.beta;
m(2,intersect(s.infectious,s.mdr)) = r.beta_mdr;


m(:,setdiff(s.infectious,[s.Ia s.Is])) = m(:,setdiff(s.infectious,[s.Ia s.Is])).*p.kappa;
m(:,intersect(s.infectious,s.hivpositive))=m(:,intersect(s.infectious,s.hivpositive)).*r.bRRhiv;

M.lambda = sparse(m);

% -- - Get the mortality rates
%HIV negative
m = zeros(1,i.nstates);
m(:) = r.mort;
m(s.hu)=r.mort_h;
m(s.ht)= r.mort_ht;
m(s.tbmort)=r.muTB;
m(s.Tx)=r.fl_mort;
m(s.Tx2)=r.sl_mort;
m(intersect(s.hu, s.tbmort))=m(intersect(s.hu, s.tbmort))*r.RRmortTBhiv;

M.mortvec = m';