
function plot_targets(runs,datatable,location)

set(0,'DefaultFigureWindowStyle','docked');

if(~isscalar(runs))
    endyr=2020;
else
    endyr=runs.endyr;
end

x=(endyr+1)-size(runs(1).popu,2):1:endyr;
% datatable{:,'mort_tbhp'} =datatable{:,'mort_tbhp'}./10;
data=datatable{location,:};
%data(5:6)=data(5:6)./10; %Readjust mort
labels=datatable.Properties.VariableNames;
limi=data([1 3 4 5 6 7 8 9]);

fields={'inc_all',	'inc_mdr',	'inc_tbhiv',	'mort_tbhn',	'mort_tbhp','notif_all'};
titles={'Incidence', 'Incidence MDR ','Incidence TB/HIV',...
    'Mortality HIV-','Mortality HIV+','Notification','HIV cascade','Notified proportions'};


ts=size(runs(1).popu,2);
nruns=1;
if (numel(runs)>1)
    runtype='mcmc';
    nruns=numel(runs);
else
    runtype='mle';
end
%functions
fun1 = @(s,field) reshape([s.(field)],ts,nruns)';

fs=11;

f1= figure;
alpha=0.2;
red = [216 23 37] ./ 255;
grey=[0.9,0.9,0.9];
lw=1.5;
xlimits=[2011 endyr];
set(f1,'color','w');

if (strcmp(location,'Ky')); coun="Kenya";
elseif (strcmp(location,'SA')); coun="South Africa";
elseif (strcmp(location,'In')); coun="India";
elseif (strcmp(location,'Br')); coun="Brazil";
end


titletxt = sprintf('%s','Model Calibration,',' ',coun);
sgt = sgtitle(titletxt,'Color','black');
sgt.FontSize = 18;

for ii=1:numel(fields)
    
    subplot(2,4,ii)
    title(titles{ii});
    hold on;
    y=fun1(runs,fields{ii});
    
    lh=plot(x,y,'linewidth',lw,'Color',grey);
    yy=prctile(y,[2.5,50,97.5],1);
    ciplot(yy(1,:),yy(3,:),x,red,alpha);
    p1=plot(x,yy(2,:),'Color',red,'linewidth',lw+0.5);
    
    iddat=find(strcmp(labels,fields{ii}));
    
    if ( (~isnan(data(iddat))) | (isempty(iddat)==1) )
        
        if strcmp(fields{ii},'inc_all')
            p1=errorbar(2012,data(iddat+1),0.35*data(iddat+1),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);
            p2=errorbar(2019,data(iddat),0.35*data(iddat),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);
        else
            p3=errorbar(2019,data(iddat),0.35*data(iddat),'o','Color','k','MarkerFaceColor','k','MarkerSize',5);
        end
        
        
    end
    set(gca,'fontsize',fs,'fontweight','bold')
    if (ii==1|| ii==4)
        ylabel('Rate per 100K','fontsize',10,'fontweight','bold')
    end
    xlim(xlimits);
    ylim([0 1.5*max([limi(ii) max(y(:))])])
    
    box on;
end



fields={'pr_onart',	'pr_onipt'};
ylabs={'%','Proportion'};
iilast=ii;
for gg=1
    ii=iilast+gg;
    subplot(2,4,ii)
    b=bar(zeros(1,2));
    names = {'onART','onIPT'};
    set(gca,'XTickLabel',[]);
    set(gca,'XTickLabel',names,'fontsize',fs);
    set(gca,'XTickLabelRotation',45)
    hold on;
    title(titles{ii});
    distr=[runs.(fields{gg,1});runs.(fields{gg,2})]';
    means=mean(distr,1);
    
    for jj=1:2
        for j=1:size(distr,1)
            
            lh=plot(jj +[-0.05 0.05],distr(j,jj).*[1 1],'linewidth',lw);
            lh.Color=[red,0.3];
        end
        p=plot(jj+[-0.05 0.05],means(jj).*[1 1],'Color','k','linewidth',lw+1);
        
        iddat=find(strcmp(labels,fields{gg,jj}));
        if ~isnan(data(iddat))
            datain=1;
            p1=errorbar(jj+0.2,data(iddat),0.2,'o','Color','k','MarkerFaceColor','k','MarkerSize',4);
            plot(jj+0.2,data(iddat),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
        end
    end
    if (datain==1)
        h=legend([p1,lh,p],'Data','Model Run','Model Mean');
        
    else
        h=legend([lh,p],'Model Run','Model Mean');
    end
    pos=get(h,'Position');
    newpos=[0.64 0.2 pos(3:4)];
    set(h,'Fontsize',fs-3,'position',newpos);
    set(gca,'fontsize',fs,'fontweight','bold')
    xlim([0 3]); ylim([0 max(distr(:))*1.5]);
    ylabel(ylabs{gg},'fontsize',10,'fontweight','bold');
    box on;
end



if strcmp(location,'Ky')
    plhv=1.6e6;
elseif strcmp(location,'SA')
    plhv=7.7e6;
elseif strcmp(location,'In')
    plhv=2.1e6;
elseif strcmp(location,'Br')
    plhv=0.9e6;
end


%% Save figure
orient(f1,'landscape')
save_fig(f1,'fits',location,runtype);

txt= sprintf('%s','results/model_fits_',coun);
print(f1,txt,'-dpdf','-r0','-bestfit')


end


