function save_results(object,type,location,runtype,folder)


filename = sprintf('%s',type,'_',location,'_',runtype,'.mat');


if (nargin<=4)
f=  fullfile('results',filename);
else
f = fullfile(folder,filename);
end
save(f,'object');

