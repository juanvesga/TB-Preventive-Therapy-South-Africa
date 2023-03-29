clear ;close all;
set(0,'DefaultFigureWindowStyle','docked')
%%

folder="results";
type  ="who";
population_output="all"; %Selectes the output as in or inc_tbhiv
hiv=strcmp(population_output,'hiv');
cuminc={sprintf('%s','cuminc','_tbhiv'*hiv)};
inc={sprintf('%s','inc','_tbhiv'*hiv)};

out=1; % 1= Cum rel reduct; 2=nnt; 3= perct reduct; 4= cases averted

% Load results
con='SA';
file=sprintf('%s',folder,'/','scenarios_',con,'_',type,'.mat');
load(file)
sa=object;


% popu;ation in 1970 for scaling up
pop70=[1.12E+07,2.28E+07,5.53E+08,5.90E+07];

sa.ired=(1-(sa.(inc{1})(:,end,2)./sa.(inc{1})(:,1,2)))*100;
sa.iav=pop70(2).*(sa.(cuminc{1})(:,end,1)-sa.(cuminc{1})(:,end,2));
arrSA=prctile(sa.arr,[2.5,50,97.5],1);
arr6HSA=prctile(sa.arr6H,[2.5,50,97.5],1);


% Labels
country={'South Africa','Kenya','India','Brazil'};
location= ["SA", "Ky" , "In" , "Br"];


outcome={'ARR','NNT','Inc. reduction(%)','Inf. Averted'};
col=linspecer(5);
roj=col(2,:); ver=col(3,:);nar=col(4,:);

%% Incidence trajectories SA
lw=1;
st=sa;
y0=prctile(st.inc(:,:,1),[2.5,50,97.5],1);
y1=prctile(st.inc(:,:,2),[2.5,50,97.5],1);
y2=prctile(st.inc(:,:,3),[2.5,50,97.5],1);
y3=prctile(st.inc(:,:,4),[2.5,50,97.5],1);

m0=prctile(st.mrt(:,:,1),[2.5,50,97.5],1);
m1=prctile(st.mrt(:,:,2),[2.5,50,97.5],1);
m2=prctile(st.mrt(:,:,3),[2.5,50,97.5],1);
m3=prctile(st.mrt(:,:,4),[2.5,50,97.5],1);

red1=round(prctile(st.arr(:,1),[2.5,50,97.5],1).*100);
red2=round(prctile(st.arr(:,2),[2.5,50,97.5],1).*100);
red3=round(prctile(st.arr(:,3),[2.5,50,97.5],1).*100);

tx1=sprintf('%s',num2str(red1(2)),'%',' ','(',...
    num2str(red1(1)),'-',num2str(red1(3)),')' );

tx2=sprintf('%s',num2str(red2(2)),'%',' ','(',...
    num2str(red2(1)),'-',num2str(red2(3)),')' );

tx3=sprintf('%s',num2str(red3(2)),'%',' ','(',...
    num2str(red3(1)),'-',num2str(red3(3)),')' );

redm1=round(prctile(st.arr_mrt(:,1),[2.5,50,97.5],1).*100);
redm2=round(prctile(st.arr_mrt(:,2),[2.5,50,97.5],1).*100);
redm3=round(prctile(st.arr_mrt(:,3),[2.5,50,97.5],1).*100);

tm1=sprintf('%s',num2str(redm1(2)),'%',' ','(',...
    num2str(redm1(1)),'-',num2str(redm1(3)),')' );

tm2=sprintf('%s',num2str(redm2(2)),'%',' ','(',...
    num2str(redm2(1)),'-',num2str(redm2(3)),')' );

tm3=sprintf('%s',num2str(redm3(2)),'%',' ','(',...
    num2str(redm3(1)),'-',num2str(redm3(3)),')' );

x=2020:2019+numel(sa.inc(1,:,1));
figure;

%%
subplot(1,2,1)
hold on
ciplot(y0(1,:),y0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,y0(2,:),'Color',col(1,:),'linewidth',lw+0.5);
box on
xlim([2020 2035])
ylim([0 max(y0(:))*1.05]);
legend([p0],{'Baseline'});
ylabel({'Rate per 100K'; 'population'});
title('TB Incidence')

%%
subplot(1,2,2)
hold on
ciplot(m0(1,:),m0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,m0(2,:),'Color',col(1,:),'linewidth',lw+0.5);
box on
xlim([2020 2035])
ylim([0 max(m0(:))*1.05]);
legend([p0],{'Baseline'});
title('TB Mortality')
%%
figure;
subplot(1,2,1)
hold on
ciplot(y0(1,:),y0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,y0(2,:),'Color',col(1,:),'linewidth',lw+0.5);

ciplot(y1(1,:),y1(3,:),x,col(2,:));alpha(0.2)
p1=plot(x,y1(2,:),'Color',col(2,:),'linewidth',lw+0.5);

box on
xlim([2020 2035])
ylim([0 max(y0(:))*1.05]);
legend([p0 p1],{'Baseline','PLHIV only'})
ylabel({'Rate per 100K'; 'population'});
title('TB Incidence')
% % annotation('textbox', [0.465, 0.98, 0.1, 0], 'string', 'Cases averted','fontsize',10,'fontweight','b','linestyle','none')
annotation('textbox', [0.465, 0.31, 0.1, 0], 'string', tx1,'fontsize',11,'linestyle','none')


%%
subplot(1,2,2)
hold on
ciplot(m0(1,:),m0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,m0(2,:),'Color',col(1,:),'linewidth',lw+0.5);

ciplot(m1(1,:),m1(3,:),x,col(2,:));alpha(0.2)
p1=plot(x,m1(2,:),'Color',col(2,:),'linewidth',lw+0.5);
title('TB Mortality')
box on
xlim([2020 2035])
ylim([0 max(m0(:))*1.05]);
legend([p0 p1],{'Baseline','PLHIV only'})
% annotation('textbox', [0.91, 0.98, 0.1, 0], 'string', 'Deaths averted','fontsize',10,'fontweight','b','linestyle','none')
annotation('textbox', [0.91, 0.31, 0.1, 0], 'string', tm1,'fontsize',11,'linestyle','none')

%%
figure;
subplot(1,2,1)
hold on
ciplot(y0(1,:),y0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,y0(2,:),'Color',col(1,:),'linewidth',lw+0.5);

ciplot(y1(1,:),y1(3,:),x,col(2,:));alpha(0.2)
p1=plot(x,y1(2,:),'Color',col(2,:),'linewidth',lw+0.5);

ciplot(y3(1,:),y3(3,:),x,col(4,:));alpha(0.2)
p3=plot(x,y3(2,:),'Color',col(4,:),'linewidth',lw+0.5);
box on
xlim([2020 2035])
ylim([0 max(y0(:))*1.05]);
legend([p0 p1 p3],{'Baseline','PLHIV only','PLHIV & HH'})
ylabel({'Rate per 100K'; 'population'});
title('TB Incidence')
% annotation('textbox', [0.465, 0.98, 0.1, 0], 'string', 'Cases averted','fontsize',10,'fontweight','b','linestyle','none')
annotation('textbox', [0.465, 0.31, 0.1, 0], 'string', tx1,'fontsize',11,'linestyle','none')
annotation('textbox', [0.465, 0.20, 0.1, 0], 'string', tx3,'fontsize',11,'linestyle','none')

%%
subplot(1,2,2)
hold on
ciplot(m0(1,:),m0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,m0(2,:),'Color',col(1,:),'linewidth',lw+0.5);

ciplot(m1(1,:),m1(3,:),x,col(2,:));alpha(0.2)
p1=plot(x,m1(2,:),'Color',col(2,:),'linewidth',lw+0.5);

ciplot(m3(1,:),m3(3,:),x,col(4,:));alpha(0.2)
p3=plot(x,m3(2,:),'Color',col(4,:),'linewidth',lw+0.5);

title('TB Mortality')
box on
xlim([2020 2035])
ylim([0 max(m0(:))*1.05]);
legend([p0 p1 p3],{'Baseline','PLHIV only','PLHIV & HH'})
% annotation('textbox', [0.91, 0.98, 0.1, 0], 'string', 'Deaths averted','fontsize',10,'fontweight','b','linestyle','none')
annotation('textbox', [0.91, 0.31, 0.1, 0], 'string', tm1,'fontsize',11,'linestyle','none')
annotation('textbox', [0.91, 0.20, 0.1, 0], 'string', tm2,'fontsize',11,'linestyle','none')



%%
% All countries
figure;
for c=1:4
if (c==1) st=sa; end


y0=prctile(st.inc(:,:,1),[2.5,50,97.5],1);
y1=prctile(st.inc(:,:,2),[2.5,50,97.5],1);
y2=prctile(st.inc(:,:,3),[2.5,50,97.5],1);
y3=prctile(st.inc(:,:,4),[2.5,50,97.5],1);

m0=prctile(st.mrt(:,:,1),[2.5,50,97.5],1);
m1=prctile(st.mrt(:,:,2),[2.5,50,97.5],1);
m2=prctile(st.mrt(:,:,3),[2.5,50,97.5],1);
m3=prctile(st.mrt(:,:,4),[2.5,50,97.5],1);

red1=round(prctile(st.arr(:,1),[2.5,50,97.5],1).*100);
red2=round(prctile(st.arr(:,2),[2.5,50,97.5],1).*100);
red3=round(prctile(st.arr(:,3),[2.5,50,97.5],1).*100);

tx1=sprintf('%s',num2str(red1(2)),'%',' ','(',...
    num2str(red1(1)),'-',num2str(red1(3)),')' );

tx2=sprintf('%s',num2str(red2(2)),'%',' ','(',...
    num2str(red2(1)),'-',num2str(red2(3)),')' );

tx3=sprintf('%s',num2str(red3(2)),'%',' ','(',...
    num2str(red3(1)),'-',num2str(red3(3)),')' );


if (c==1)

% annotation('textbox', [0.465, 0.98, 0.1, 0], 'string', 'Cases averted','fontsize',10,'fontweight','b','linestyle','none')
elseif (c==2)
% annotation('textbox', [0.91, 0.98, 0.1, 0], 'string', 'Cases averted','fontsize',10,'fontweight','b','linestyle','none')
    
end


x=2020:2019+numel(sa.inc(1,:,1));

subplot(2,2,c)
hold on
ciplot(y0(1,:),y0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,y0(2,:),'Color',col(1,:),'linewidth',lw+0.5);

% ciplot(y1(1,:),y1(3,:),x,col(2,:));alpha(0.2)
% p1=plot(x,y1(2,:),'Color',col(2,:),'linewidth',lw+0.5);

% ciplot(y2(1,:),y2(3,:),x,col(3,:));alpha(0.2)
% p2=plot(x,y2(2,:),'Color',col(3,:),'linewidth',lw+0.5);

ciplot(y3(1,:),y3(3,:),x,col(4,:));alpha(0.2)
p3=plot(x,y3(2,:),'Color',col(4,:),'linewidth',lw+0.5);

xlim([2020 2035])
ylim([0 max(y0(:))*1.1])

if (c==4)
legend([p0, p3],{'Baseline','PLHIV & 100% HH'})
end
if (c==1||c==3)
ylabel({'Incidence Rate'; 'per 100K population'});
end
title(country{c})
box on
if (c==1)
% annotation('textbox', [0.465, 0.73, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.465, 0.68, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.465, 0.65, 0.1, 0], 'string', tx3,'fontsize',11,'linestyle','none')
elseif(c==2)
% annotation('textbox', [0.91, 0.68, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.91, 0.64, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.91, 0.61, 0.1, 0], 'string', tx3,'fontsize',11,'linestyle','none')

elseif(c==3)
% annotation('textbox', [0.465, 0.37, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.465, 0.26, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.465, 0.22, 0.1, 0], 'string', tx3,'fontsize',11,'linestyle','none')

elseif(c==4)
% annotation('textbox', [0.91, 0.34, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.91, 0.31, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.91, 0.27, 0.1, 0], 'string', tx3,'fontsize',11,'linestyle','none')

end


end


%%
% All countries
h=figure;
fs=12;
for c=1
if (c==1) st=sa; end


y0=prctile(st.inc(:,:,1),[2.5,50,97.5],1);
y1=prctile(st.inc(:,:,2),[2.5,50,97.5],1);
y2=prctile(st.inc(:,:,3),[2.5,50,97.5],1);
y3=prctile(st.inc(:,:,4),[2.5,50,97.5],1);

m0=prctile(st.mrt(:,:,1),[2.5,50,97.5],1);
m1=prctile(st.mrt(:,:,2),[2.5,50,97.5],1);
m2=prctile(st.mrt(:,:,3),[2.5,50,97.5],1);
m3=prctile(st.mrt(:,:,4),[2.5,50,97.5],1);

red1=round(prctile(st.arr(:,1),[2.5,50,97.5],1).*100);
red2=round(prctile(st.arr(:,2),[2.5,50,97.5],1).*100);
red3=round(prctile(st.arr(:,3),[2.5,50,97.5],1).*100);

tx1=sprintf('%s',num2str(red1(2)),'%',' ','(',...
    num2str(red1(1)),'-',num2str(red1(3)),')' );

tx2=sprintf('%s',num2str(red2(2)),'%',' ','(',...
    num2str(red2(1)),'-',num2str(red2(3)),')' );

tx3=sprintf('%s',num2str(red3(2)),'%',' ','(',...
    num2str(red3(1)),'-',num2str(red3(3)),')' );


if (c==1)

% annotation('textbox', [0.465, 0.98, 0.1, 0], 'string', 'Cases averted','fontsize',10,'fontweight','b','linestyle','none')
elseif (c==2)
% annotation('textbox', [0.91, 0.98, 0.1, 0], 'string', 'Cases averted','fontsize',10,'fontweight','b','linestyle','none')
    
end


x=2020:2019+numel(sa.inc(1,:,1));


hold on
ciplot(y0(1,:),y0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,y0(2,:),'Color',col(1,:),'linewidth',lw+0.5);

 ciplot(y1(1,:),y1(3,:),x,col(2,:));alpha(0.2)
 p1=plot(x,y1(2,:),'Color',col(2,:),'linewidth',lw+0.5);

% ciplot(y2(1,:),y2(3,:),x,col(3,:));alpha(0.2)
% p2=plot(x,y2(2,:),'Color',col(3,:),'linewidth',lw+0.5);

ciplot(y3(1,:),y3(3,:),x,col(4,:));alpha(0.2)
p3=plot(x,y3(2,:),'Color',col(4,:),'linewidth',lw+0.5);

xlim([2020 2035])
ylim([0 max(y0(:))*1.1])

if (c==4)
legend([p0, p1, p3],{'Baseline','PLHIV only','PLHIV & 100% HH'}, 'location',...
    'SE','FontSize',fs)
end
if (c==1||c==3)
ylabel({'Incidence rate'; 'per 100K population'});
end
title(country{c})
box on
if (c==1)
 annotation('textbox', [0.465, 0.76, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.465, 0.68, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.465, 0.69, 0.1, 0], 'string', tx3,'fontsize',10,'linestyle','none')
elseif(c==2)
 annotation('textbox', [0.91, 0.77, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.91, 0.64, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.91, 0.7, 0.1, 0], 'string', tx3,'fontsize',10,'linestyle','none')

elseif(c==3)
 annotation('textbox', [0.465, 0.39, 0.1, 0], 'string', "<1%",'fontsize',10,'linestyle','none')
% annotation('textbox', [0.465, 0.26, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.465, 0.28, 0.1, 0], 'string', tx3,'fontsize',10,'linestyle','none')
xlabel("Year")
elseif(c==4)
 annotation('textbox', [0.91, 0.38, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.91, 0.31, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.91, 0.35, 0.1, 0], 'string', tx3,'fontsize',10,'linestyle','none')
xlabel("Year")
end


end


% 
% 
% set(h,'PaperOrientation','landscape');
% 
% txt= sprintf('%s','C:\Users\JFV09\Dropbox\Shared Folders\TBshared_Nim_Juan\Future PT paper\Figures\fig2');
% print(h,txt,'-dpdf','-r0','-bestfit')
% 

%%
% check 6H
h=figure;
fs=12;
for c=1
if (c==1) st=sa; end
if (c==2) st=ky; end
if (c==3) st=in; end
if (c==4) st=br; end

y0=prctile(st.inc(:,:,1),[2.5,50,97.5],1);
y6=prctile(st.inc(:,:,6),[2.5,50,97.5],1);


m0=prctile(st.mrt(:,:,1),[2.5,50,97.5],1);
m6=prctile(st.mrt(:,:,6),[2.5,50,97.5],1);


red6=round(prctile(st.arr(:,5),[2.5,50,97.5],1).*100);


tx1=sprintf('%s',num2str(red6(2)),'%',' ','(',...
    num2str(red6(1)),'-',num2str(red6(3)),')' );


if (c==1)

% annotation('textbox', [0.465, 0.98, 0.1, 0], 'string', 'Cases averted','fontsize',10,'fontweight','b','linestyle','none')
elseif (c==2)
% annotation('textbox', [0.91, 0.98, 0.1, 0], 'string', 'Cases averted','fontsize',10,'fontweight','b','linestyle','none')
    
end


x=2020:2019+numel(sa.inc(1,:,1));


hold on
ciplot(y0(1,:),y0(3,:),x,col(1,:));alpha(0.2)
p0=plot(x,y0(2,:),'Color',col(1,:),'linewidth',lw+0.5);

 ciplot(y6(1,:),y6(3,:),x,col(2,:));alpha(0.2)
 p1=plot(x,y6(2,:),'Color',col(2,:),'linewidth',lw+0.5);

xlim([2020 2035])
ylim([0 max(y0(:))*1.1])

if (c==4)
legend([p0, p1, p3],{'Baseline','Expand 6H'}, 'location',...
    'SE','FontSize',fs)
end
if (c==1||c==3)
ylabel({'Incidence rate'; 'per 100K population'});
end
title(country{c})
box on
if (c==1)
 annotation('textbox', [0.465, 0.76, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.465, 0.68, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.465, 0.69, 0.1, 0], 'string', tx3,'fontsize',10,'linestyle','none')
elseif(c==2)
 annotation('textbox', [0.91, 0.77, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.91, 0.64, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.91, 0.7, 0.1, 0], 'string', tx3,'fontsize',10,'linestyle','none')

elseif(c==3)
 annotation('textbox', [0.465, 0.39, 0.1, 0], 'string', "<1%",'fontsize',10,'linestyle','none')
% annotation('textbox', [0.465, 0.26, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.465, 0.28, 0.1, 0], 'string', tx3,'fontsize',10,'linestyle','none')
xlabel("Year")
elseif(c==4)
 annotation('textbox', [0.91, 0.38, 0.1, 0], 'string', tx1,'fontsize',10,'linestyle','none')
% annotation('textbox', [0.91, 0.31, 0.1, 0], 'string', tx2,'fontsize',8,'linestyle','none')
annotation('textbox', [0.91, 0.35, 0.1, 0], 'string', tx3,'fontsize',10,'linestyle','none')
xlabel("Year")
end


end
