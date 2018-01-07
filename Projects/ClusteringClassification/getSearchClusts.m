clear all; close all;

reload = 0;
reClust = 1;
doClip = 0;
areaCrit = {'FEF','F2','SEF','SC'};
recCrit = {'rich','brad','plex','tdt'};
blWind = -300:-100;
doReward = 0;
minN = 10;
saveLoad = 1;

matLoc = [mlRoot,'Dropbox/Schall-Lab/dataMats/search'];

if reload
    [searchSDF,searchTimes,searchAreas,searchRecSys, searchRow, searchMonks, searchOtherTasks, searchSess, searchChans] = klPullAllSDFs('search','-a',areaCrit,'-r',doReward,'rec',recCrit,'-c',doClip);
    if saveLoad
        thisDate = datestr(today,25);
        save([matLoc,filesep,'searchSDFs-',thisDate(~ismember(thisDate,'/')),'-r',num2str(doReward),'c',num2str(doClip),'.mat'],'search*');
    end
else
%     load('searchSDFs-171015.mat');
    % Get most recent searchSDFs file from matLoc
    matFiles = dir([matLoc,'/*r',num2str(doReward),'c',num2str(doClip),'.mat']);
    if isempty(matFiles)
        [searchSDF,searchTimes,searchAreas,searchRecSys, searchRow, searchMonks, searchOtherTasks, searchSess, searchChans] = klPullAllSDFs('search','-a',areaCrit,'-r',doReward,'rec',recCrit,'-c',doClip);
        if saveLoad
            thisDate = datestr(today,25);
            save([matLoc,filesep,'searchSDFs-',thisDate(~ismember(thisDate,'/')),'-r',num2str(doReward),'c',num2str(doClip),'.mat'],'search*');
        end
    else
        matDates = cellfun(@(x) str2double(x(12:17)),{matFiles.name});
        load([matLoc,filesep,matFiles(matDates==max(matDates)).name]);
    end
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
myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preVis,visTrans,visSust,preMov,postMov,nextVis};
myEpocInds = [1,1,1,2,2,2,3,3,3,4,4,4];
myEpocWinds = {[-200:300],[-300:200],[-200:300],[-300:200]};
if doReward
    myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preRwd,earlyRwd,lateRwd,preVis,visTrans,visSust,preMov,postMov,nextVis,preRwd,earlyRwd,lateRwd};
    myEpocInds = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6];
    myEpocWinds = {[-200:300],[-300:200],[-100:100],[-200:300],[-300:200],[-100:100]};%myEpocWinds = cat(2,myEpocWinds,{[-100:100],[-100:100]});
end

inAreas = searchAreas(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit));

% Convert those SDFs to a cell that can be input to the clustering
inSDFs{1} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),1),'UniformOutput',0));
inSDFs{2} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),2),'UniformOutput',0));
% inSDFs{3} = cell2mat(cellfun(@(x) x(1,:)-x(3,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),1),'UniformOutput',0));
% inSDFs{4} = cell2mat(cellfun(@(x) x(1,:)-x(3,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),2),'UniformOutput',0));
inSDFs{3} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),1),'UniformOutput',0));
inSDFs{4} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),2),'UniformOutput',0));
normSD = [1:4];
if doReward
    inSDFs{5} = inSDFs{4};
    inSDFs{4} = inSDFs{3};
    inSDFs{3} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),3),'UniformOutput',0));
%     inSDFs{6} = cell2mat(cellfun(@(x) x(1,:)-x(3,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),3),'UniformOutput',0));
    inSDFs{6} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),3),'UniformOutput',0));
    normSD = cat(2,normSD,[5,6]);
end

plotSDFs{1} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),1),'UniformOutput',0));
plotSDFs{2} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),2),'UniformOutput',0));
plotSDFs{3} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),1),'UniformOutput',0));
plotSDFs{4} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),2),'UniformOutput',0));
if doReward
    plotSDFs{5} = plotSDFs{4};
    plotSDFs{4} = plotSDFs{3};
    plotSDFs{3} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),3),'UniformOutput',0));
    plotSDFs{6} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),3),'UniformOutput',0));
end

% Get time cell
inTimes{1} = nanmean(cell2mat(searchTimes(:,1)),1);
inTimes{2} = nanmean(cell2mat(searchTimes(:,2)),1);
inTimes{3} = inTimes{1};
inTimes{4} = inTimes{2};
if doReward
    inTimes{5} = inTimes{4};
    inTimes{4} = inTimes{3};
    inTimes{3} = nanmean(cell2mat(searchTimes(:,3)),1);
    inTimes{6} = inTimes{3};
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
inAreas = inAreas(goodRows);

if reClust
    [sortIDs, idxDist, raw, respSumStruct, respLink] = klDoMetaClusteringv2(inSDFs,inTimes,'-e',myEpocs,'-ei',myEpocInds,'-er',myEpocWinds,'-r',0,'sd',normSD,'-mn',minN);
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

normResp = klNormRespv2(plotSDFs,inTimes,'ztr','-r',myEpocWinds,'bl',blWind,'sd',1:4);%normSD);

plotSearchClusts(normResp,inTimes,sortIDs,respAdj,raw,myK);



