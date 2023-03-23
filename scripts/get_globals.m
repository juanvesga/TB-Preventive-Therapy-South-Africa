

% Get all necessary model specifications, fixed parameters data
% and function handles

location = "SA";

% Some running parameters
forg= 0.5; %
mdr_effect= 0.5;
reinf_protect=1;
hiv_switch=1; %turns hiv on (1) or off(0)
correct_hiv_eff = 1;% 0.7;% 0.78;
post_size=1000;
endyear=2035;
fityear=2025; %2022 To fit HLM targets, 2023 for WHO
scenario     = "results"; % res or full (full is extreme values);res_split_pot
baseline     = "normal";     %Alternative is 6H
RifFree      =  0 ;% Default=0; 1 is Rifampin free regimen (no DR BArrier)
drmin        = 0.95; % Minimal boundary of the DR acquisition
forg         = 0.5;%frg; % Level of completion to which forgiveness apply( default is 50%)
pot_split    = 1;% Potency only after some level of complition (def is 1, no splfit)
fup          = 2;    %Follow-up time for efficacy cohort (defaut is 2)
attr         = 7;  % Set to 6 for No_Rif scenario (def is 7)
mdr_effect   = 0.5; %effect of regimen on rif strains (def is 50%)(vary 25-75)
pt_type      = "who";
range        = 'optimal';
llk          = 'llk'; % llk or lsq
ptpop        = "all"; %plhv simulates interventions only among HIV on ART
pt_scale     = 2;% years of pt clae up (default=2, vary 10 and 15)
%% Check and override conditions according to selected sceanrios in master script
if (exist('ptpop','var'))
    if(strcmp(ptpop,"plhv"))
        scenario=sprintf('%s','hiv',scenario);
    end
else
    scenario='sens';
    ptpop='all';
end





% Options: "IPT", "3HP", "1HP" <-- Not used like this now
regimen="6IPT";


% Load Data
data  =load_data(strcat(pwd,'/data/WHO_data.xlsx'));

T=readtable('/data/WHO_data_ranges.xlsx');
id=find(strcmp(T{:,1},location));
ranges=[T{id,4}./T{id,3} T{id,3}./T{id,3} T{id,5}./T{id,3} ];
tmp=reshape(T{:,3},size(data,2),size(data,1));
for kk=1:size(data,2)
    data.(kk)=tmp(kk,:)';
end

pars  =load_data('data/programatic_parameters.xlsx');
hivrts=load_data('data/UNAIDS_hiv_rates.xlsx');
hlm_targ  =load_data('data/HLM_tpt_targets.xlsx');
tmp       =load_data('data/WHO_tpt_targets.xlsx');
who_targ=tmp(:,1);
reg_profile=load_data('data/tpt_regimen_profiles.xlsx');
house_comp  =load_data('data/country_household_size.xlsx');
hivpoints=hivrts{location,:};

%% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Map model structure in a system of indexes for constructing model
% equations
gps.strain = {'ds' ,'mdr'};
gps.hiv    = {'hn' ,'hu','ht'};
states0 = {'U','Pu1','Pu2','Qu','Ru','Qc','Rc'};
states1 = {'Lf','Ls','Pf1','Pf2','Ps1','Ps2','Qf','Qs','Rf','Rs','Ia','Is','E','Dx','Tx','Tx2','Rlo','Rhi','R'};

[i, s, d, lim] = get_addresses({states0, gps.hiv}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states1,gps.strain, gps.hiv }, i, s, d, lim);
d = char(d);

%Include the auxiliaries (these create model outputs)
auxnames = {'inc', 'notif' , 'mort' ,'pt' , 'newart','ptpmo','daly','ptcomp','ptforg','remote','recent'};
auxinds  = [  5     , 3    ,  5     ,  3  ,    2    ,   2    ,  4  ,    2,       2,     1      ,   1];
for ii = 1:length(auxnames)
    inds = lim + (1:auxinds(ii));
    i.aux.(auxnames{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;

% dictionary of model states and useful aggregations
s.nstates=(1:i.nstates);
s.hivpositive=[s.hu, s.ht];
s.infectious = [s.Ia,s.Is,s.E,s.Dx,intersect(s.Tx,s.mdr)];
s.tbmort     = [s.Ia,s.Is,s.E,s.Dx,s.Tx, s.Tx2];
s.prevalent  = unique([s.infectious, s.Tx,s.Tx2]);
s.tbcare     = [s.Dx s.Tx s.Tx2];
s.notbcare   = [s.Pu1, s.Pf1, s.Ps1, s.Pu2,s.Pf2, s.Ps2 , ...
    s.Qu, s.Qc, s.Qf, s.Qs,...
    s.Ru,s.Rc,s.Rf,s.Rs,s.Ia,s.Is,s.E, s.Rhi, s.Rlo s.R];
s.TBmort     = s.infectious;
s.ipt        = [s.Pu1 s.Pf1 s.Ps1];
s.ipt_all     = [s.Pu1 s.Pf1 s.Ps1 s.Pu2 s.Pf2 s.Ps2];
s.sought= [s.Dx s.E ];
s.treated= [s.Tx, s.Tx2];
s.postdx =[s.Dx s.Tx s.Tx2 s.Rlo s.Rhi s.R];
s.ptcomp =[s.Qu, s.Qc, s.Qf, s.Qs ];


% Parameters to be calubrated with ranges
xi = [];

xnames =...
    {'beta',...             % Transmission rate per capita (per year)
    'beta_mdr',...          % Same but for DR
    'bRRhiv',...            % RR of TB progression given in PLHIV
    'reinfhiv',...          % Reinfection factor in PLHIV
    'fast',...              % yearly pc rate of fast progression to TB
    'symp_del',...          % yearly rate of developing symptoms for asymp TB
    'careseek',...          % Care seeking rate per year
    'ARTrec',....           % ART recruitment rate
    'IPThiv',...            % IPT coverage among those starting ART
    'Dxnhiv',...            % TB Diagnosis probability in non HIV
    'Dxhiv',...             % TB Diagnosis probability in PLHIV
    'Txinit',...            % Probability of TB treatment initiation once diagnosed
    'muTB',...              % TB yearly mortality rate
    'RRmortTBhiv',...       % RR of TB mortality in those PLHIV
    'ntpcov'};              % Overall coverage of the TB programme
xnums  = ones(1,size(xnames,2));
lim = 0;
for ii = 1:length(xnames)
    inds = lim + (1:xnums(ii));
    xi.(xnames{ii}) = inds;
    lim = inds(end);
end
bds=zeros(length(xnames),2);
xi.nx = lim;
bds(xi.beta,:)             = [0.1 50];
bds(xi.beta_mdr,:)         = [0.1 50];
bds(xi.bRRhiv,:)           = [0.6 1];
bds(xi.reinfhiv,:)         = [1 12];
bds(xi.fast,:)             = [0.01 1.2];
bds(xi.symp_del,:)         = [0.1 15];
bds(xi.careseek,:)         = [0.5 12];
bds(xi.ARTrec,:)           = [0.1 10];
bds(xi.IPThiv,:)           = [0.01 1];
bds(xi.Dxnhiv,:)           = [0.5 0.98];
bds(xi.Dxhiv,:)            = [0.5 0.98];
bds(xi.Txinit,:)           = [0.5 0.98];
bds(xi.muTB,:)             = [0.08 0.2];
bds(xi.RRmortTBhiv,:)      = [1 6];
bds(xi.ntpcov,:)           = [0.1 1];




prm.bds = bds';

%-- Call distributions for lhd
% data_llk=data;
% data_llk.mort_tbhp=data_llk.mort_tbhp.*10;
% if (strcmp(llk,'llk'))
%
%     lhd=Make_distribution_fns(data_llk,ranges,location,'epi');
%     lhd_hlm=Make_distribution_fns(hlm_targ,[],location,'hlm');
%     lhd_who=Make_distribution_fns(who_targ,[],location,'who');
% else
%     dat= data_llk;% data{location,:};
%     dat_hlm=hlm_targ{location,:};
%     lhd=    @(x)  -sum( (dat-x).^2);
%     lhd_hlm= @(x)  -( (dat_hlm-x).^2);
% end

% Selectors for the computing incidence (useful for selecting/aggregating
% during matrix multiplication
% This avoids aggregating transitions not allowed in model
tmp = ones(i.nstates);
tmp(s.hn,s.hu)=0;
tmp(s.hu,s.hn)=0;
tmp(s.hn,s.ht)=0;
tmp(s.ht,s.hn)=0;
tmp(s.hu,s.ht)=0;
tmp(s.ht,s.hu)=0;
check= sparse(tmp - diag(diag(tmp)));
sel.check=check;

%  Incidence: 1.Total, 2.MDR, 3.TB/HIV,  4. <15 , 5.>=15 , 6. MDR_acq , HIV Incidence
tmp = zeros(4,i.nstates);
tmp(1,s.Ia) = 1;                           %All
tmp(2,intersect(s.Ia,s.mdr)) = 1;          %MDR
tmp(3,intersect(s.Ia,s.hivpositive)) = 1;%HIV/TB
tmp(4,intersect(intersect(s.Ia,s.hivpositive),s.mdr)) = 1; %HIV/MDR
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Ia,:) = 1;
tmp=tmp.*check;
sel.inc = sparse(tmp - diag(diag(tmp)));

% Aggregate for incidence cohort
tmp = zeros(2,i.nstates);
tmp(1,intersect(s.Ia,[s.hn,s.hu])) = 1;      %All
tmp(2,intersect(s.Ia,s.ht)) = 1;%HIV/TB
agg.inc_c = sparse(tmp);

tmp = zeros(i.nstates);
tmp(intersect(s.Tx,s.mdr),intersect(s.Tx,s.ds)) = 1;
tmp=tmp.*check;
sel.acqu = sparse(tmp - diag(diag(tmp)));

tmp = zeros(i.nstates);
tmp([s.Lf, s.Qf, s.Pf1, s.Pf2, s.Rf],[s.U s.Ls s.Rlo s.Rhi s.R s.Pu1 s.Pu2 s.Ps1 s.Ps2 s.Qu s.Qc s.Qs s.Ru s.Rc s.Rs]) = 1;
tmp=tmp.*check;
sel.ltbi = sparse(tmp - diag(diag(tmp)));

% TB notification
tmp = zeros(3,i.nstates);
tmp(1,s.treated)  = 1;
tmp(2,s.Tx2)  = 1;
tmp(3,intersect(s.treated,s.hivpositive))  = 1;

agg.notif = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.treated,s.Dx) = 1;
tmp=tmp.*check;
sel.notif = sparse(tmp - diag(diag(tmp)));


% Selectors mortality
tmp = zeros(5,i.nstates);
tmp(1,s.nstates) = 1;
tmp(2,s.tbmort) = 1;
tmp(3,intersect(s.tbmort,s.hivpositive)) = 1;
tmp(4,intersect(s.tbmort,s.mdr)) = 1;
tmp(5,intersect(intersect(s.tbmort,s.hivpositive),s.mdr)) = 1;
agg.mort = sparse(tmp);


% IPT
tmp = zeros(3,i.nstates);
tmp(1,s.ipt)  = 1;
tmp(2,intersect(s.ipt,s.ht))  = 1;% PT
tmp(3,intersect(s.ipt,[s.hn s.hu]))  = 1;% PT a
agg.ipt = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.ipt,[s.U s.Lf s.Ls]) = 1;
tmp=tmp.*check;
tmp(intersect(s.ipt,s.ht),intersect(s.hu,[s.U s.Lf s.Ls])) = 1;
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
sel.ipt = sparse(tmp - diag(diag(tmp)));

% New on ART
tmp = zeros(2,i.nstates);
tmp(1,[intersect(s.ht,[s.U s.Ls s.Lf]) intersect(s.ht,s.ipt)]  ) = 1;
tmp(2,s.ht)  = 1;
agg.art = sparse(tmp);

tmp = zeros(i.nstates);
tmp=tmp.*check;
tmp(s.ht,s.hu) = 1;
sel.art = sparse(tmp - diag(diag(tmp)));


%PT PMO
tmp = zeros(2,i.nstates);
tmp(1,s.ipt_all)  = 1;
tmp(2,intersect(s.ipt_all,s.ht))  = 1;% PT amongst HIV on care
agg.ptpmo = sparse(tmp);


%PT PMO
tmp = zeros(2,i.nstates);
tmp(1,s.infectious)  = 1;
tmp(2,intersect(s.infectious,s.hivpositive))  = 1;% PT amongst HIV on care
agg.daly = sparse(tmp);


% IPT completion and forgiveness
tmp = zeros(2,i.nstates);
tmp(1,s.ptcomp)  = 1;
tmp(2,intersect(s.ptcomp,s.ht))  = 1;% PT
agg.ptcomp = sparse(tmp);

tmp = zeros(i.nstates);
tmp(intersect(s.ptcomp,s.mdr),intersect(s.ds,[ s.Pf2, s.Ps2])) = 1;
tmp=tmp.*check;
tmp(s.ptcomp,[s.Pu2, s.Pf2, s.Ps2]) = 1;
tmp([s.Ia],[s.Pf1, s.Ps1,s.Pf2, s.Ps2]) = 1;
tmp([s.Qc],[s.Pf1, s.Ps1,s.Pf2, s.Ps2]) = 1;
tmp(s.ds,s.mdr)=0;
tmp(s.mdr,s.ds)=0;
tmp(intersect(s.ptcomp,s.ht),intersect(s.hu,[s.Pu2, s.Pf2, s.Ps2])) = 1;
sel.ptcomp = sparse(tmp - diag(diag(tmp)));


% Remote Incidence
tmp = zeros(1,i.nstates);
tmp(1,s.Ia) = 1;                           %All
agg.remote = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Ia,[s.Ls,s.Ps1,s.Ps2,s.Qs,s.Rs]) = 1;
tmp=tmp.*check;
sel.remote = sparse(tmp - diag(diag(tmp)));

% Recent Incidence
tmp = zeros(1,i.nstates);
tmp(1,s.Ia) = 1;                           %All
agg.recent = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Ia,[s.Lf,s.Pf1,s.Pf2,s.Qf,s.Rf]) = 1;
tmp=tmp.*check;
sel.recent = sparse(tmp - diag(diag(tmp)));

%% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Model parameters


% Linear rates of transition
r.reactivation  =  0.001;
r.reactivation_fast =  1/2;
r.careseeking2  = [0 0 0];
r.access        = 0;

% Diagnosis stage
r.Dx = 52;
p.Dx = [0,0,0];%by hiv status
p.Tx_init = 0.75;
p.smear_sens=0.8; %Smear test speccificity (Swai F 2011, BMC resreach notes)
p.smear_spec=0.94;%Smear test speccificity (Swai F 2011, BMC resreach notes)
p.xpert_sens=0.9; %Smear test speccificity (Swai F 2011, BMC resreach notes)
p.xpert_spec=0.99;%Smear test speccificity (Swai F 2011, BMC resreach notes)
p.smear= 0.38;%prob of microscopy by sector
p.xpert= 0.17;%prob of microscopy by sector
p.yrntp = pars{location,'yrntp'};% care seeking 0.18pry as estimated by Olney et al 2016
p.yrart = pars{location,'yrhaart'};
p.ntpcov=0.5; % Scale NTP in phases..50%
r.cs2           =  2;
r.csRRhiv       = 3;


% Treatment stage (Usign country specific WHO database for Tx outcomes)
r.Tx = 2;
p.fl_fail = pars{location,'fl_fail'};
p.fl_lost = pars{location,'fl_lost'};
p.fl_die  = pars{location,'fl_die'};
p.fl_suc  =  1-(p.fl_fail+p.fl_lost+p.fl_die);%

r.TxOut   = r.Tx/(p.fl_suc+p.fl_fail); % total rate of treatment outcomes
r.default = r.TxOut*p.fl_lost;
r.fl_mort = r.TxOut*p.fl_die;
p.cure    = p.fl_suc/(p.fl_suc+p.fl_fail);

%MDR TB
r.MDR_acqu = 0.05;
p.SL_trans = 0.8; %Transfer from FL
r.Tx2      = 0.5;
p.sl_fail = pars{location,'sl_fail'};
p.sl_lost = pars{location,'sl_lost'};
p.sl_die  = pars{location,'sl_die'};
p.sl_suc  = 1-(p.sl_fail+p.sl_lost+p.sl_die);%pars{location,'sl_succ'};

r.TxOut2   = r.Tx2/(p.sl_suc+p.sl_fail); % total rate of treatment outcomes
r.default2 = r.TxOut2*p.sl_lost;
r.sl_mort = r.TxOut2*p.sl_die;
p.cure2    = p.sl_suc/(p.sl_suc+p.sl_fail);


% Natural history parameters
p.kappa     = 1;
r.selfcure      =  0.16;
r.muTB          =  0.16;
r.RRmortTBhiv   = 3;
r.relapse   = [0.032 0.14 0.0015];
r.fast_react      =0.0826;%Menxies [2.4090 ,  0.9855 , 0.0985 , 0.0985 ];      % Ragonnet 2018_Epidemics
r.slow_react      =0.000594;%Menzies [6.9350e-09 ,   0.0023 , 0.0012 , 0.0012 ]; % Ragonnet 2018_Epidemics
r.slow            =0.8728;%Menzies   [4.38, 4.38 , 1.9710 , 1.9710];             % Ragonnet 2018_Epidemics


p.imm       = 0.5;
p.crossg    = 0.3;


% HIV cascades
r.progRRhiv = 26; % Getahun , Selwyn
p.hivtest_notb = pars{location,'hivtest_notb'};% care seeking 0.18pry as estimated by Olney et al 2016
p.hivtest_tb   = pars{location,'hivtest_tb'}  ;% WHO report % known status.. care seeking inside TB care system
p.tbtest_hiv   = pars{location,'tbtest_hiv'}  ;% Tb test among HIV pos
p.art_notb     = pars{location,'art_notb'}   ;% Coverage expected in 2015 * proportion estimated to adhreing and getting suppressed
p.art_tb       = pars{location,'art_tb'}   ;% Coverage expected in TB+ 2015 * proportion estimated to adhreing and getting suppressed
r.art_dropout  = pars{location,'art_dropout'};% Average UNGAS2016 + AMPATH/
r.hivdecline   = pars{location,'hivdecline'};    % Annual Decline in HIV incidence after 2017
p.xpert        = pars{location,'xpert'};
r.reinfhiv=1; % Susceptibility to reinfection on HIV+
r.ARTred    = 0.35;        % Reduction in progression given ART (Suthar2010PlosMed)
r.ARTrec    = 0;           % Recruitment probability
p.hivcq_ad      = 1;


%Preventive therapy
p.IPT=zeros(1,3) ;%ia,ig
p.forg=forg;
r.ptdur           =  reg_profile{regimen,'ptdur'};    % Duration of protection post regimen
r.ptregdur        =  reg_profile{regimen,'ptregdur'};     % Regimen duration (halved to show first and 2nd half)
p.pteffi          =  reg_profile{regimen,'pteffi'};     % efficacy
p.ptdrbarrier     =  reg_profile{regimen,'ptdrbarr'};     % DR Barrier
p.ptforg          =  reg_profile{regimen,'ptforg'};     % Regimen Forgiveness
p.ptcomp          =  reg_profile{regimen,'ptcomp'};     % Regimen completion (ease of adherence
p.nomido          =  reg_profile{regimen,'ptnomido'};    %
p.midoforg        =  reg_profile{regimen,'ptmidoforg'};    %
p.potency         =  reg_profile{regimen,'potency'};    %
p.reinfprotection =  reinf_protect;
r.corr_hiv  = correct_hiv_eff;%  0.75; % ratio of efficacy HIV+/HIV-

r.regdur1= r.ptregdur*p.forg;
r.regdur2= r.ptregdur*(1-p.forg);
r.pt_resist_acq =  (12./[r.regdur1 r.regdur2])*(1-p.ptdrbarrier)./p.ptdrbarrier; %Rate of acquisition of resistance during PT
r.ptdefault     =  (12/(r.regdur1+ r.regdur2))*(1-p.ptcomp)/p.ptcomp; % estimated time lost of protection due to completion rates
r.outReg        =  12./([r.regdur1 r.regdur2]);
r.outIPT        =  12/(r.ptdur);      % 24 months protection Rangaka-Martins etc and 3HP
p.effondr       =  mdr_effect; % Effect of PT on Dr strains

p.housedist =[1 1 1]; % rate control for household distr of TST
p.LTBItest  = 0;% set to 1 to set on PT only latent infections
% Demographic terms
p.lex              = pars{location,'lex'};% Life expectancy
r.mort             = 1/p.lex;%pars{location,{'mrt1','mrt2','mrt3','mrt4'}};   % Non-disease-related mortality from WHO Life Tables nMx
r.mort_h = 1/15;          % Non-TB-related mortality in untreated HIV+
r.mort_ht= 1/36;          % Non-TB-related mortality in treated HIV+ (assumption: ex6tebeds live as far as observed, 2017-1982)
r.mort_TB = 1/6;          % Untreated(10%)/pre-treatment(90%)TB death in HIV negative
r.mort_TBtx = 1/6;        % Untreated(10%)/pre-treatment(90%)TB death in HIV negative
r.mort_TBh= 1/2;          % Untreated TBH with HIV+
r.mort_TBh_tx =1/2;       % Death during treatment in HIV+
r.RRmortTBhiv =1.8;
p.growth            = pars{location,'growth'};
p.frac_pop          = pars{location,{'fpop1','fpop2','fpop3','fpop4'}};        % Population fraction <15 2016;
r.ageing =[ 1/3.5, 1/7, 1/90];          % Rate of ageing into >15 (tweaked to match simulated population fraction)
p.pop1970=pars{location,'pop1970'};
p.house_age = house_comp{location, 8:end };

% Interventions
r.OptimART = 1;           % Parameter to optimize HIV cascade in intervention
r.acf_asym=0;
r.acf_sym= 0;
r.nutri=0;

% DALYs
p.yll=[68.7-10 68.7-30]; % Life Years lost in chidrenand aults
p.weights=[0.33 1];      % Utility weights: 1. Untreated TB 2.Deceased


% Switch Model Control
p.priv_off   = 1; %(==0 cuts pass to private)
r.turnoffHIV = hiv_switch; %(==0 cuts HIV)
p.ARToff     = 1; %(==0 turns ARToff)
r.cohort     =1; % 0=on, cuts FOI
r.lam        =zeros(2,1);%foi for cohort
r.pot_split=pot_split;
r.fup=fup; %Efficacy measured at x yeras of follow-up

prm.p = p; prm.r = r; prm.endy=endyear;prm.fity=fityear; prm.ptscale=pt_scale;
ref.i = i; ref.s = s; ref.d = d; ref.xi = xi;
if exist('drmin','var')
    ref.strdr='perm';
    if drmin>0.95
        ref.strdr=sprintf('%s','perm',string(drmin*100));
    end
else
    drmin=0.95;
end

prm.ptbase=[9,0,0.95,0.1,0.617];
prm.switchhiv=strcmp(ptpop,"plhv");
% Create instances of model function for different return types
obj     = @(x)           get_objective(x, prm, ref, sel, agg, gps, lhd,hivpoints,endyr);
obj_mcmc = @(x)  getfield(get_objective(x, prm, ref, sel, agg, gps, lhd,hivpoints,endyr),'llk');
%     obj_spx = @(x)  -getfield(get_objective(x, prm, ref, sel, agg, gps, lhd,hivpoints),'llk');

% --- Parameterization of the model


%LHS PT params boundaries
ptbds=zeros(10,2);
ptbds(1,:)=[0 120];   % Durabilty
ptbds(2,:)=[1 3];     % Duration (mo)
ptbds(3,:)=[0 1];     % Suppression
ptbds(4,:)=[drmin 1]; % MDR acquit
ptbds(5,:)=[0.5 0.8]; % Forgiveness
ptbds(6,:)=[0.8 0.9]; % Ease of adherence
ptbds(7,:)=[0 .99];   % Potency
ptbds(8,:)=[0 120];   % Durabilty 6H
ptbds(9,:)=[0 1]; % Suppression 6H
ptbds(10,:)=[0 0.99]; % Potency 6H

if (strcmp(scenario,"res_full"))
    
    ptbds= ptbds.*[0.75 1.25];
    ptbds(3,2)=1;
    ptbds(4,2)=1;
    ptbds(5,2)=1;
    ptbds(6,2)=1;
    ptbds(7,2)=1;
    ptbds(8,:)=[0 120]; % Durabilty 6H
    ptbds(9,:)=[0 1]; % Suppression 6H
    ptbds(10,:)=[0 0.99]; % Potency 6H
    
    
elseif(strcmp(scenario,"res_lowcompletion"))
    
    ptbds=zeros(10,2);
    ptbds(1,:)=[0 120];   % Durabilty
    ptbds(2,:)=[1 3];     % Duration (mo)
    ptbds(3,:)=[0 1];     % Suppression
    ptbds(4,:)=[drmin 1]; % MDR acquit
    ptbds(5,:)=[0.5 0.8]; % Forgiveness
    ptbds(6,:)=[0.4 0.42]; % Ease of adherence
    ptbds(7,:)=[0 .99];   % Potency
    ptbds(8,:)=[0 120];   % Durabilty 6H
    ptbds(9,:)=[0 1]; % Suppression 6H
    ptbds(10,:)=[0 0.99]; % Potency 6H
    
end



