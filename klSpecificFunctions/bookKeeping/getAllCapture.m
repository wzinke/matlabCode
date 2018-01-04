monk        = {'Gauss','Helmholtz'};
xlFile      = './cosDataBookKeeping_capture.xlsx';

for im = 1:2
    
    [excelNum,excelText,excelAll] = xlsread(xlFile,monk{im});
    excelNum = cat(1,nan(1,size(excelNum,2)),excelNum);
    excelNum = cat(2,nan(size(excelNum,1),size(excelAll,2)-size(excelNum,2)),excelNum);
    hRow     = find(strcmpi(excelAll(:,1),'Session'));
    
    if im == 1, 
        allMonkChans = (hRow+1):size(excelNum,1);
    else
        allMonkChans = (hRow+1):size(excelNum,1);
    end
    close all; 
    klPlotCapturev4(allMonkChans,'-smoc',1,monk{im},0,1);
    
    close all;
end
    