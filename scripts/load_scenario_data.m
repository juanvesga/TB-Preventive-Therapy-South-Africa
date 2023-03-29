%% This script loads necessary data for running a scenario analysis 

% File paths 
f1=sprintf('%s',"results",'/','output','_',location,'_','mcmc','.mat');
f2=sprintf('%s',"results",'/','house_pars','_',location,'_','mcmc','.mat');
f3=sprintf('%s',"results",'/','whofits','_',location,'_','mcmc','.mat');
f4=sprintf('%s',"results",'/','lhs_',location,'_who.mat');

% Load calibrated model simulatiomn results
load(f1);  
x   =object.x;

% Get latest state of the model
for qq=1:post_size  
    sfin1(qq,:)=object(qq).sfin;
end
sfin=sfin1;              
id= find(strcmp(xnames,'IPThiv'));

% Load household composition parameters 
load(f2);
prm.hcpars=object.pars;
prm.hcdist=object.dist;

%  Load parameters for baseline fit to WHO TPT targets
load(f3);
prm.id= find(strcmp(xnames,'IPThiv'));
ptpars=object.xpt;

% Load analysis results of TPT regimens 
load(f4)
runs=object;