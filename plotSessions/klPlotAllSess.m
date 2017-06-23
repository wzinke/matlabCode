% klPlotAllSess
close all;

monk = {'Gauss','Helmholtz'};
saveSess = 1;

for im = 1:length(monk),

    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});

    uSess = unique(excelAll(strcmp(cellfun(@class,excelAll(:,2),'UniformOutput',0),'char'),2));

    uSess(ismember(uSess,{'Name','Session'})) = [];

    for is = 1:length(uSess),
        klPlotSessv1(uSess{is},'-m',monk{im},'-s',saveSess);
        if saveSess,
            close all;
        end
    end
end