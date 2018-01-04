clear all; close all;

reload = 0;
reClust = 1;
areaCrit = {'SEF'};
recCrit = {'rich','brad','plex','tdt'};
blWind = -300:-100;
myN = 10;
metaN = 10;

if reload
    [sefSDF,sefTimes,sefAreas,sefRecSys] = klPullAllSDFs('mg','-a',areaCrit,'-r',1);
else
    load('sefMGSDFs-171015.mat');
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
myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preRwd,earlyRwd,lateRwd};
myEpocInds = [1,1,1,2,2,2,3,3,3];
myEpocWinds = {[-200:300],[-300:200],[-100:100]};

% Convert those SDFs to a cell that can be input to the clustering
inSDFs{1} = sefSDF{1}(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),:);
inSDFs{2} = sefSDF{2}(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),:);
inSDFs{3} = sefSDF{3}(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),:);

plotSDFs{1} = sefSDF{1}(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),:);
plotSDFs{2} = sefSDF{2}(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),:);
plotSDFs{3} = sefSDF{3}(ismember(sefAreas,areaCrit) & ismember(sefRecSys,recCrit),:);

% Get time cell
inTimes{1} = nanmean(cell2mat(sefTimes(:,1)),1);
inTimes{2} = nanmean(cell2mat(sefTimes(:,2)),1);
inTimes{3} = nanmean(cell2mat(sefTimes(:,3)),1);

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
    [sortIDs, idxDist, raw, respSumStruct, respLink] = klDoMetaClusteringv2(inSDFs,inTimes,'-e',myEpocs,'-ei',myEpocInds,'-er',myEpocWinds,'-r',0,'-n',myN,'-mn',metaN);
end
respAdj = respLink; respAdj(:,3) = respLink(:,3)-min(respLink(:,3));
for ic = 1:size(raw,2)
    for ir =  ic:size(raw,1)
        raw(ir,ic) = raw(ic,ir);
    end
end

nBad = sum(isnan(sortIDs(:,1)));
nGood = sum(~isnan(sortIDs(:,1)));

for ik = 1:size(sortIDs,2)
    nNan(ik) = (sum(isnan(sortIDs(:,ik)))-nBad)./nGood;
end

myK = find(nNan <= .1,1,'last');

normResp = klNormRespv2(plotSDFs,inTimes,'ztr','-r',myEpocWinds,'bl',blWind);

plotMGClusts(normResp,inTimes,sortIDs,respAdj,raw,myK);



