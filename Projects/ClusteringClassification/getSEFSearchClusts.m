clear all; close all; 

reload = 0;
reClust = 1;
areaCrit = {'SEF'};
recCrit = {'rich','brad'};%,'plex','tdt'};
blWind = -300:-100;
doReward = 0;
if reload
    [sefSDF,sefTimes,sefAreas,sefRecSys] = klPullAllSDFs('search','-a',areaCrit,'-r',1);
else
    load('sefSearchSDFs-171019.mat');
end

% Set analysis epochs
preVis = -100:0;
visTrans = 50:100;
visSust = 100:150;
preMov = -50:0;
postMov = 0:50;
nextVis = 50:100;
preRwd = -50:0;
earlyRwd = 0:50;
lateRwd = 50:100;
myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preRwd,earlyRwd,lateRwd,preVis,visTrans,visSust,preMov,postMov,nextVis,preRwd,earlyRwd,lateRwd};
myEpocInds = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6];
myEpocWinds = {[-200:300],[-300:200],[-200:100],[-200:300],[-300:200],[-200:100]};

normSD = [1,2,3,4,5,6];

% Convert those SDFs to a cell that can be input to the clustering
inSDFs{1} = cell2mat(cellfun(@(x) x(1,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),1),'UniformOutput',0));
inSDFs{2} = cell2mat(cellfun(@(x) x(1,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),2),'UniformOutput',0));
inSDFs{3} = cell2mat(cellfun(@(x) x(1,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),3),'UniformOutput',0));
% inSDFs{4} = cell2mat(cellfun(@(x) x(1,:)-x(3,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),1),'UniformOutput',0));
% inSDFs{5} = cell2mat(cellfun(@(x) x(1,:)-x(3,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),2),'UniformOutput',0));
% inSDFs{6} = cell2mat(cellfun(@(x) x(1,:)-x(2,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),3),'UniformOutput',0));
inSDFs{4} = cell2mat(cellfun(@(x) x(3,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),1),'UniformOutput',0));
inSDFs{5} = cell2mat(cellfun(@(x) x(3,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),2),'UniformOutput',0));
inSDFs{6} = cell2mat(cellfun(@(x) x(2,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),3),'UniformOutput',0));


plotSDFs{1} = cell2mat(cellfun(@(x) x(1,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),1),'UniformOutput',0));
plotSDFs{2} = cell2mat(cellfun(@(x) x(1,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),2),'UniformOutput',0));
plotSDFs{3} = cell2mat(cellfun(@(x) x(1,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),3),'UniformOutput',0));
plotSDFs{4} = cell2mat(cellfun(@(x) x(3,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),1),'UniformOutput',0));
plotSDFs{5} = cell2mat(cellfun(@(x) x(3,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),2),'UniformOutput',0));
plotSDFs{6} = cell2mat(cellfun(@(x) x(2,:),sefSDF(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),3),'UniformOutput',0));

% Get time cell
inTimes{1} = nanmean(cell2mat(sefTimes(:,1)),1);
inTimes{2} = nanmean(cell2mat(sefTimes(:,2)),1);
inTimes{3} = nanmean(cell2mat(sefTimes(:,3)),1);
inTimes{4} = inTimes{1};
inTimes{5} = inTimes{2};
inTimes{6} = inTimes{3};

if ~doReward
    inSDFs = inSDFs([1,2,4,5]);
    plotSDFs = plotSDFs([1,2,4,5]);
    inTimes = inTimes([1,2,4,5]);
    myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preVis,visTrans,visSust,preMov,postMov,nextVis};
    myEpocInds = [1,1,1,2,2,2,3,3,3,4,4,4];
    myEpocWinds = {[-200:300],[-300:200],[-200:300],[-300:200]};
    normSDF = 1:4;
end
    

catNorms = [];
for i = 1:length(inSDFs)
    catNorms = cat(2,catNorms,inSDFs{i}(:,ismember(inTimes{i},myEpocWinds{i})));
end
goodRows = sum(isfinite(catNorms),2) == size(catNorms,2);

for i = 1:length(inSDFs)
    inSDFs{i} = inSDFs{i}(goodRows,:);
    plotSDFs{i} = plotSDFs{i}(goodRows,:);
end

if reClust
    [sortIDs, idxDist, raw, respSumStruct, respLink] = klDoMetaClusteringv2(inSDFs,inTimes,'-e',myEpocs,'-ei',myEpocInds,'-er',myEpocWinds,'-r',0,'sd',normSD);
end
respAdj = respLink; respAdj(:,3) = respLink(:,3)-min(respLink(:,3));
for ic = 1:size(raw,2)
    for ir =  ic:size(raw,1)
        raw(ir,ic) = raw(ic,ir);
    end
end

nBad = sum(isnan(sortIDs(:,1)));
nGood = sum(~isnan(sortIDs(:,1)));
nNan = nan(1,size(sortIDs,2));
for ik = 1:size(sortIDs,2)
    nNan(ik) = (sum(isnan(sortIDs(:,ik)))-nBad)./nGood;
end

myK = find(nNan <= .1,1,'last');

normResp = klNormRespv2(plotSDFs,inTimes,'ztr','-r',myEpocWinds,'bl',blWind);

plotSearchClusts(normResp,inTimes,sortIDs,respAdj,raw,myK);



