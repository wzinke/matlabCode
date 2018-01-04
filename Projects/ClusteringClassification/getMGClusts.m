clear all; close all;

reload = 1;
reClust = 1;
areaCrit = {'FEF','F2','SEF','SC'};%{'FEF'};%,'F2'};
recCrit = {'rich','brad','plex','tdt'};
blWind = -300:-100;
myN = 10;
metaN = 10;
doReward = 0;
saveLoad = 1;
doClip = 1;

matLoc = '~/Dropbox/Schall-Lab/dataMats/memGuide';

if reload
    [mgSDF,mgTimes,mgAreas,mgRecSys,mgRows,mgMonks,mgOtherTasks,mgSess,mgChans,mgParams] = klPullAllSDFs('mg','-r',doReward,'rec',recCrit,'-a',areaCrit,'-c',doClip);
    if saveLoad
        thisDate = datestr(today,25);
        save([matLoc,filesep,'mgSDFs-',thisDate(~ismember(thisDate,'/')),'-r',num2str(doReward),'.mat'],'mg*');
    end
else
%     load('mgSDFs-171025.mat');
%     load('mgSDFs-171015.mat');
    % Get most recent mgSDFs file from matLoc
    matFiles = dir([matLoc,'/*r',num2str(doReward),'.mat']);
    matDates = cellfun(@(x) str2double(x(8:13)),{matFiles.name});
    load([matLoc,filesep,matFiles(matDates==max(matDates)).name]);
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
myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis};
myEpocInds = [1,1,1,2,2,2];
myEpocWinds = {[-200:300],[-300:200]};
if doReward
    myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preRwd,earlyRwd,lateRwd};
    myEpocInds = [1,1,1,2,2,2,3,3,3];
    myEpocWinds = {[-200:300],[-300:200],[-200:200]};
end
inAreas = mgAreas(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit));
inRecSys = mgRecSys(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit));

% Convert those SDFs to a cell that can be input to the clustering
inSDFs{1} = mgSDF{1}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
inSDFs{2} = mgSDF{2}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
if doReward
    inSDFs{3} = mgSDF{3}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
end

plotSDFs{1} = mgSDF{1}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
plotSDFs{2} = mgSDF{2}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
if doReward
    plotSDFs{3} = mgSDF{3}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
end

% Get time cell
inTimes{1} = nanmean(cell2mat(mgTimes(:,1)),1);
inTimes{2} = nanmean(cell2mat(mgTimes(:,2)),1);
if doReward
    inTimes{3} = nanmean(cell2mat(mgTimes(:,3)),1);
end

catNorms = [];
for i = 1:length(inSDFs)
    catNorms = cat(2,catNorms,inSDFs{i}(:,ismember(inTimes{i},myEpocWinds{i})));
end
goodRows = sum(isfinite(catNorms),2) == size(catNorms,2);

inAreas = inAreas(goodRows);

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



