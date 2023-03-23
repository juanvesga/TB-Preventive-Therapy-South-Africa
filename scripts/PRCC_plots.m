
set(0,'DefaultFigureWindowStyle','docked')
%%
scen=scenario;
type  ="who";
baseline= "";
func = @(params,outcome) prcc(params,outcome);
titxt = {'Partial Rank Correlation (PRCC) for TB cases averted'};
analysispop="all";%selectes the interventions as only among HIV on ARt or all
population_output="all"; %Selectes the output as in or inc_tbhiv
hiv=strcmp(population_output,'hiv');
cuminc={sprintf('%s','cuminc','_tbhiv'*hiv)};
inc={sprintf('%s','inc','_tbhiv'*hiv)};
out=4; % 1= Cum rel reduct; 2=nnt; 3= perct reduct; 4= cases averted
gridset={'grid2'};
nboots=200; % Number of bootstsraps

ptbds(1,:)=[0 1200];% Durabilty (mo)
ptbds(2,:)=[1 3];% Duration (mo)
ptbds(3,:)=[0.7 1];% Efficacy
ptbds(4,:)=[0.95 1]; % MDR acquit
ptbds(5,:)=[0.5 0.8]; % Forgiveness
ptbds(6,:)=[0.8 0.9]; % Ease of adherence
ptbds(7,:)=[0 1]; % potency
if (strcmp(scenario,"res_lowcompletion"))
    ptbds(6,:)=[0.4 0.42];
end

% Load results
con='SA';
file=sprintf('%s',scen,'_','lhs',baseline,'_',con,'_who.mat');
load(file)
sa=object;

pop70=[1.12E+07,2.28E+07,5.53E+08,5.90E+07];

sa.ired=(1-(sa.(inc{1})(:,end)./sa.(inc{1})(:,1)))*100;
sa.iav=pop70(2).*(sa.(cuminc{1})(1,end)-sa.(cuminc{1})(:,end));
sa.iavs=pop70(2).*(sa.(cuminc{1})(1,:)-sa.(cuminc{1})(:,:));
S=[sa.(gridset{1}),sa.arr(2:end),sa.nnt(2:end),sa.ired(2:end),sa.iav(2:end)];
All=S;
G1=sa.grid;
effes=sa.effs(2:end,end);
effarts=sa.effarts(2:end,end);

% Labels
country={'South Africa','Kenya','India','Brazil'};
location= [ "SA" ,"Ky", "In" , "Br"];
names={'Durability','Regimen duration', 'Efficacy','DR-Barrier','Forgiveness','Ease-of-adherence','Cure'};
names2={'Regimen duration', 'Efficacy','DR-Barrier','Forgiveness','Ease-of-adherence'};
names3={'Regimen duration', 'Efficacy','DR-Barrier','Forgiveness','Ease-of-adherence','Cure','Waning','Suppression'};
names4={'Regimen duration','DR-Barrier','Forgiveness','Ease-of-adherence','Cure','Waning','Suppression'};
names5={'Proportion cure','Durability of non-curative protection','Strength of non-curative protection'};
names6={'Regimen duration','DR-Barrier','Forgiveness','Ease-of-adherence','Cure','Suppression'};


outcome={'ARR','NNT','Inc. reduction(%)','Inf. Averted'};
col=linspecer(5);
roj=col(2,:); ver=col(3,:);nar=col(4,:);
natt=numel(names);


%% Plot PRCCs
c=1;
f1=figure;
set(gcf,'color','w');
pars_id=[2 3 4 5 6];
labs=names(pars_id);
natt_ou=numel(labs);
nco=1;
ij=0;
[med, lo, up]=deal(zeros(nco,natt_ou));
params=squeeze(All(:,pars_id));
z1=squeeze(All(:,[1 7]));
z2=squeeze(G1(:,[1 3]));
z=[z1,z2];
iq=find(params(:,2)>ptbds(3,1));
params=params(iq,:);
z=z(iq,:);
outcome=squeeze(All(iq,natt+out));
[bootstat,bootsam] = bootstrp(nboots,func,params,outcome);
yy=prctile(bootstat,50);
yy_low=prctile(bootstat,2.5);
yy_up=prctile(bootstat,97.5);
med(:)=yy;
lo(:)=yy_low;
up(:)=yy_up;
[~, ij]=sort(yy);
i=ij;

hold on;
b=barh(yy(i),'FaceAlpha',1,'EdgeColor','none');
er = errorbar(yy(i),1:natt_ou,yy(i)-yy_low(i),yy_up(i)-yy(i),'.','horizontal');
er.LineWidth = 1.5;
er.Color = 'k';
er.MarkerSize = 1;
ylim([0.5 natt_ou+.5])
xlim([0 1]);
yticks([1 2 3 4 5])
yticklabels(labs(i))
set(b, 'FaceColor',col(c,:));
xlabel('PRCC');
t = title(country{c}, 'Units', 'normalized', 'Position', [0.7,0.15,1]);
box on
grid
set(gca,'fontsize',12);


%% Scatter plots of epidemic impact vs regimen attributes
f2=figure;
pars_id=[2 3 4 5 6];
labs=names(pars_id);
natt_ou=numel(labs);
nco=1;
[med, lo, up]=deal(zeros(nco,natt_ou));

tiledlayout('flow')
params=squeeze(All(:,pars_id));
z1=squeeze(All(:,[1 7]));
z2=squeeze(G1(:,[1 3]));
z=[z1,z2];
iq=find(params(:,2)>ptbds(3,1));
params=params(iq,:);
z=z(iq,:);
outcome=squeeze(All(iq,natt+out));

[bootstat,bootsam] = bootstrp(nboots,func,params,outcome);
yy=prctile(bootstat,50);
yy_low=prctile(bootstat,2.5);
yy_up=prctile(bootstat,97.5);

med(c,:)=yy;
lo(c,:)=yy_low;
up(c,:)=yy_up;
[~, ij]=sort(yy);
i=ij;

nm=labs(flip(i));


for g=1:natt_ou
    nexttile;
    ii=flip(i);
    x=rescale(params(:,ii(g)));
    y=rescale(outcome);
    plot(x,y,'.','Color',col(c,:)*1,'markersize',2);
    
    ylim([min(0,min(y)) max(y)])
    set(gca,'fontsize',8);
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    ax=gca;
    ax.YRuler.Exponent = 0;
    xticklabels([]);
    lag=0.8;
    title(nm{g});
    lag=1;
    xlabel({'Property Value';' (normalized)'});
    ylabel({'Epidemic';'impact'});
    yl=get(gca,'ylim');
    xl=get(gca,'xlim');
        
    % Annotate PRCC
    yl=get(gca,'ylim');
    xl=get(gca,'xlim');
    text(xl(1)+0.2,yl(2)-0.3,{'PRCC=',num2str(prccs(c,ii(g),2))}, 'fontsize',7);
    
    box on
end



