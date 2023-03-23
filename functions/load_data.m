function out=load_data(address)



[nums, txt, ~] =xlsread(address);

rowNames = txt(2:end,1);

if (size(txt,2)>1)
    
    colNames = txt(1,2:size(nums,2)+1);
else
    colNames=    matlab.lang.makeValidName(cellstr(string(nums(1,:))));
    nums=nums(2:end,:);
end


sTable = array2table(nums,'RowNames',rowNames,'VariableNames',colNames);


out=sTable;