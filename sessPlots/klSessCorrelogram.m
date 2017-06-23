% klSessCorrelogram

close all;
monk = {'Gauss','Helmholtz'};
xlFile = 'klDataBookKeeping_mg.xlsx';


colStr = 'type';
areaCrit = {'FEF'};
randType = 'inSess'; % inSess, shuffle, specVals (more to come)
visualize = 1;

randReps = 1000;

allSessHist = [];
allSessCount = [];
monkRandHist = {};
monkRandCount = {};
monkID = [];
allSessX = [];
for im = 1:length(monk),
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    depthCol  = find(strcmp(excelAll(4,:),'depth (chan#)'));
    testCol = find(strcmp(excelAll(4,:),colStr));
    areaCol = find(strcmp(excelAll(4,:),'area'));
    
    if ischar(excelAll{5,testCol}),
        uTest = unique(excelAll(5:end,testCol));
    else
        uTest = unique(cell2mat(excelAll(5:end,testCol)));
        uTest(isnan(uTest)) = [];
    end
    tstCount = nan(1,length(uTest));
    for it = 1:length(uTest),
        if ischar(excelAll{5,testCol}),
            tstCount(it) = sum(ismember(excelAll(5:end,testCol),uTest{it}));
        else
            tstCount(it) = sum(ismember(excelNum(5:end,testCol),uTest(it)));
        end
    end
    
    uSess = unique(excelAll(5:end,2));
    monkSessHist = []; monkSessCount = []; monkX = [];
    nSess = 0;
    for is = 1:length(uSess),
        myChanInds = strcmp(excelAll(:,2),uSess{is});
        if sum(ismember(excelAll(myChanInds,areaCol),areaCrit)) ~= sum(myChanInds),
           
            continue;
        end
        nSess = nSess+1;
        myDepths = excelNum(myChanInds,depthCol);
        myTypes  = excelAll(myChanInds,testCol);
        
        [tmpHist, tmpCount, tmpX, randHist{nSess}, randCount{nSess}] = klDepthCorrelogram(myTypes,myDepths,'rtype',randType,tstCount,'randreps',randReps);
        monkSessHist = cat(1,monkSessHist,tmpHist);
        monkSessCount = cat(1,monkSessCount,tmpCount);
        monkX = cat(1,monkX,tmpX);
    end
    monkID = cat(1,monkID,ones(length(uSess),1).*im);
    allSessHist = cat(1,allSessHist,monkSessHist);
    allSessCount = cat(1, allSessCount, monkSessCount);
    allSessX    = cat(1,allSessX,monkX);
    
    % Restructure randHist, randCount
    newRandHist = cell(1,randReps);
    newRandCount = cell(1,randReps);
    for ir = 1:randReps,
        for is = 1:nSess,
            if ~isempty(randHist{is}),
                newRandHist{1,ir} = cat(1,newRandHist{1,ir},randHist{is}{ir});
                newRandCount{1,ir} = cat(1,newRandCount{1,ir},randCount{is}{ir});
            end
        end
    end
    
    monkRandHist = cat(1,monkRandHist,newRandHist);
    monkRandCount = cat(1,monkRandCount,newRandCount);
end

allRandHist = cell(1,randReps);
allRandCount = cell(1,randReps);
endRandHist = []; endRandCount = [];
for ir = 1:randReps,
    for im = 1:size(monkRandHist,1),
        allRandHist{1,ir} = cat(1,allRandHist{1,ir},monkRandHist{im,ir});
        allRandCount{1,ir} = cat(1,allRandCount{1,ir},monkRandCount{im,ir});
    end
    endRandHist(ir,:) = sum(allRandHist{1,ir},1);
    endRandCount(ir,:) = sum(allRandCount{1,ir},1);
end

if visualize,
    figure();
    plot(nanmean(allSessX,1),nansum(allSessHist)./nansum(allSessCount),'k','linewidth',2);
    hold on;
    plot(nanmean(allSessX,1), nanmean(endRandHist./endRandCount,1), 'color', [.8 .2 .2],'linewidth',2);
    plot(nanmean(allSessX,1), nanmean(endRandHist./endRandCount,1)+nanstd(endRandHist./endRandCount,[],1),'--', 'color', [.8 .2 .2]);
    legend({'Observed','Randomized Mean','Random Mean +/- 1*SD'});
    plot(nanmean(allSessX,1), nanmean(endRandHist./endRandCount,1)-nanstd(endRandHist./endRandCount,[],1),'--', 'color', [.8 .2 .2]);
    xlabel('Channel Separation','fontsize',14);
    ylabel('% Same Category','fontsize',14);
    title(sprintf('Channel Cross-Correlation: %s',colStr));
    
    figure();
    plot(nanmean(allSessX,1),((nansum(allSessHist,1)./nansum(allSessCount,1))-nanmean(endRandHist./endRandCount,1))./nanstd(endRandHist./endRandCount,[],1),'linewidth',2,'color','k');
    hl(1) = hline(-2); hl(2) = hline(2); set(hl,'color',[.8 .2 .2],'linewidth',1,'linestyle','--')
    legend({'Z-Scored Observation'});
    xlabel('Channel Separation','fontsize',14);
    ylabel(sprintf('Z-Same Category: %d Reps',randReps),'fontsize',14);
    title(sprintf('Z-Scored Cross-Correlation: %s',colStr),'fontsize',18);
    
end
