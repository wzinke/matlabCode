%% pltGroupChans

% Some stuff here to get the excel row indices you want...
% excelRows = ...
[excelNum,excelText,excelAll] = xlsread(xlFile,'Gauss');
excelNum = cat(1,nan(1,size(excelNum,2)),excelNum);
excelNum = cat(2,nan(size(excelNum,1),size(excelAll,2)-size(excelNum,2)),excelNum);
gPosF2 = find(excelNum(:,posCol) == 1 & strcmp(excelAll(:,4),'F2'));

excelRows = gPosF2;

[allMeans, allStds, allTimes, allEvents] = klPlotChansv3(excelRows,'-mos','Gauss',0,0);

minTimes = cellfun(@min,allTimes);
minTime = min(minTimes,[],1);

numNanFront = minTimes - (repmat(minTime,size(allMeans,1),1));
[adjMeans, adjStds] = deal(cell(size(numNanFront,1),2));

for iu = 1:size(allMeans,1),
    adjMeans{iu,1} = cat(2,nan(1,numNanFront(iu,1)),allMeans{iu,1});
    adjMeans{iu,2} = cat(2,nan(1,numNanFront(iu,2)),allMeans{iu,2});
    adjStds{iu,1}  = cat(2,nan(1,numNanFront(iu,1)),allStds{iu,1});
    adjStds{iu,2}  = cat(2,nan(1,numNanFront(iu,2)),allStds{iu,2});
end

vectLengths = cellfun(@length,adjMeans);
maxLength = max(vectLengths,[],1);
numNanBack = repmat(maxLength,size(adjMeans,1),1) - vectLengths;

for iu = 1:size(allMeans,1)
    adjMeans{iu,1} = cat(2,adjMeans{iu,1},nan(1,numNanBack(iu,1)));
    adjMeans{iu,2} = cat(2,adjMeans{iu,2},nan(1,numNanBack(iu,2)));
    adjStds{iu,1}  = cat(2,adjStds{iu,1},nan(1,numNanBack(iu,1)));
    adjStds{iu,2}  = cat(2,adjStds{iu,2},nan(1,numNanBack(iu,1)));
end

grpMean{1} = nanmean(cell2mat(adjMeans(:,1)),1);
grpStd{1}  = nanstd(cell2mat(adjMeans(:,1)),1)./sqrt(size(adjMeans,1));
grpMean{2} = nanmean(cell2mat(adjMeans(:,2)));
grpStd{2}  = nanstd(cell2mat(adjMeans(:,2)),1)./sqrt(size(adjMeans,1));

figure();
pltMeanStd(1:length(grpMean{1}),grpMean{1},grpStd{1},'k');
