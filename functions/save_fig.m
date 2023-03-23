function save_fig(object,type,location,runtype)


fig = get(groot,'CurrentFigure');
if ~isempty(fig)
filename = sprintf('%s',type,'_',location,'_',runtype,'.fig');
f = fullfile('results',filename);
savefig(object,f,'compact'); 
end
