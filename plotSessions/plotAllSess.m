% plotAllSess

monks = {'Helmholtz','Gauss'};
xlFile          = './cosDataBookKeeping_mg.xlsx';

area = {'All'};

for im = 1:2
    [excelNum,excelText,excelAll] = xlsread(xlFile,monks{im});
    excelNum = cat(1,nan(1,size(excelNum,2)),excelNum);
    excelNum = cat(2,nan(size(excelNum,1),size(excelAll,2)-size(excelNum,2)),excelNum);
    
    if ismember('All',area),
        areaRows = 2:size(excelAll,1);
    else
        areaRows = find(ismember(excelAll(:,strcmpi(excelAll(1,:),'Area')),area));
    end
    allSess = unique(excelAll(areaRows,1));
    allSess(ismember(allSess,'Session')) = [];
    
    for is = 1:length(allSess)
        pltSessChans(monks{im},allSess{is});
        close all;
    end
    
end
