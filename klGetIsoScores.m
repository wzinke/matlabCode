monk = {'Gauss','Helmholtz','Darwin'};
xlFile = './klDataBookKeeping_mg.xlsx';

load masterSorts
load allIsoScores

global excelNum excelAll

doFresh = 0;

% allIsoScores = cell(1,length(monk));
for im = 1:length(monk),
    
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    
    if doFresh || isempty(allIsoScores{im}),
        startRow = 5;
        allIsoScores{im} = nan(size(excelNum,1),4);
    else
        startTmp = find(~isnan(allIsoScores{im}(:,1)),1,'last')+1;
        startRow = max([5,startTmp]);
    end
    fprintf('Monkey %d: ',im);
    if startRow > length(masterSorts.(monk{im})),
        fprintf('Already done');
    end
    for ir = startRow:length(masterSorts.(monk{im})),
        fprintf('Row %d of %d...',ir,length(masterSorts.(monk{im})));
        
        % Get substructure
        thisStruct = masterSorts.(monk{im})(ir);
        
        % Load file for waveforms
        [path,file] = klRowToFile(ir,'-m',monk{im},'-r',0,'-w',1);
        load([path,file{1}]);
        
        [snr,isoScore,fnScore,fpScore] = klGetAllIso(wave.waves,thisStruct);
        
        allIsoScores{im}(ir,:) = [snr,isoScore,fnScore,fpScore];
        
        if mod(ir,20) == 0,
%             allIsoScores{im} = isos;
            save allIsoScores.mat allIsoScores
%             fprintf('Monkey %d: Row %d of %d\n',im,ir,length(masterSorts.(monk{im})));
        end
        if ir < length(masterSorts.(monk{im})),
            for ib = 1:length(sprintf('Row %d of %d...',ir,length(masterSorts.(monk{im})))),
                fprintf('\b');
            end
        end
    end
    fprintf('\n');
    
%    allIsoScores{im} = isos;
    save allIsoScores.mat allIsoScores
end
save allIsoScores.mat allIsoScores