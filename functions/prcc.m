function out=prcc(pars,outcome,z)


if nargin==2
corrs=partialcorri(outcome,pars,'type','Spearman');
else
corrs=partialcorri(outcome,pars,z,'type','Spearman');
end    
    
    
% out=corrs(end,1:end-1);
out=corrs';
end