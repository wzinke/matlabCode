% clear all; 
close all;

reload = 0;
reClust = 1;
areaCrit = {'FEF'};
blWind = -300:-100;
doReward = 0;
myN = 10;
if reload
    [paSDF,paTimes,paTypes,paRTs] = klPullProAnti;
else
    load('paSDFs-171026.mat');
end
paAreas = cell(size(paSDF,1),1);
[paAreas{1:size(paSDF,1)}] = deal('FEF');

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
myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preVis,visTrans,visSust,preMov,postMov,nextVis,preVis,visTrans,visSust,preMov,postMov,nextVis,preVis,visTrans,visSust,preMov,postMov,nextVis};
myEpocInds = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8];
myEpocWinds = {[-200:300],[-300:200],[-200:300],[-300:200],[-200:300],[-300:200],[-200:300],[-300:200]};
if doReward
    myEpocs = cat(2,myEpocs,{preRwd,earlyRwd,lateRwd,preRwd,earlyRwd,lateRwd});
    myEpocInds = cat(2,myEpocInds,[5,5,5,6,6,6]);
    myEpocWinds = cat(2,myEpocWinds,{[-100:100],[-100:100]});
end


% Convert those SDFs to a cell that can be input to the clustering
inSDFs{1} = cell2mat(cellfun(@(x) x(1,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
inSDFs{2} = cell2mat(cellfun(@(x) x(1,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
inSDFs{3} = cell2mat(cellfun(@(x) x(3,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
inSDFs{4} = cell2mat(cellfun(@(x) x(3,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
% inSDFs{5} = cell2mat(cellfun(@(x) x(1,:)-x(2,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
% inSDFs{6} = cell2mat(cellfun(@(x) x(1,:)-x(2,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
% inSDFs{7} = cell2mat(cellfun(@(x) x(3,:)-x(4,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
% inSDFs{8} = cell2mat(cellfun(@(x) x(3,:)-x(4,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
inSDFs{5} = cell2mat(cellfun(@(x) x(2,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
inSDFs{6} = cell2mat(cellfun(@(x) x(2,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
inSDFs{7} = cell2mat(cellfun(@(x) x(4,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
inSDFs{8} = cell2mat(cellfun(@(x) x(4,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
normSD = [1:8];
if doReward
    inSDFs{5} = cell2mat(cellfun(@(x) x(1,:),paSDF(ismember(paAreas,areaCrit),3),'UniformOutput',0));
    inSDFs{6} = cell2mat(cellfun(@(x) x(1,:)-x(3,:),paSDF(ismember(paAreas,areaCrit),3),'UniformOutput',0));
    normSD = cat(2,normSD,5);
end

plotSDFs{1} = cell2mat(cellfun(@(x) x(1,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
plotSDFs{2} = cell2mat(cellfun(@(x) x(1,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
plotSDFs{3} = cell2mat(cellfun(@(x) x(2,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
plotSDFs{4} = cell2mat(cellfun(@(x) x(2,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
plotSDFs{5} = cell2mat(cellfun(@(x) x(3,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
plotSDFs{6} = cell2mat(cellfun(@(x) x(3,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
plotSDFs{7} = cell2mat(cellfun(@(x) x(4,:),paSDF(ismember(paAreas,areaCrit),1),'UniformOutput',0));
plotSDFs{8} = cell2mat(cellfun(@(x) x(4,:),paSDF(ismember(paAreas,areaCrit),2),'UniformOutput',0));
if doReward
    inSDFs{5} = cell2mat(cellfun(@(x) x(1,:),paSDF(ismember(paAreas,areaCrit),3),'UniformOutput',0));
    inSDFs{6} = cell2mat(cellfun(@(x) x(3,:),paSDF(ismember(paAreas,areaCrit),3),'UniformOutput',0));
end

% Get time cell
inTimes{1} = nanmean(cell2mat(paTimes(:,1)),1);
inTimes{2} = nanmean(cell2mat(paTimes(:,2)),1);
inTimes{3} = inTimes{1};
inTimes{4} = inTimes{2};
inTimes{5} = inTimes{1};
inTimes{6} = inTimes{2};
inTimes{7} = inTimes{1};
inTimes{8} = inTimes{2};
if doReward
    inTimes{5} = nanmean(cell2mat(paTimes(:,3)),1);
    inTimes{6} = inTimes{5};
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
    [sortIDs, idxDist, raw, respSumStruct, respLink] = klDoMetaClusteringv2(inSDFs,inTimes,'-e',myEpocs,'-ei',myEpocInds,'-er',myEpocWinds,'-r',0,'sd',normSD,'do',1,'-mn',myN);
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

normResp = klNormRespv2(plotSDFs,inTimes,'ztr','-r',myEpocWinds,'bl',blWind,'sd',1:8);%normSD);

plotPAClusts(normResp,inTimes,sortIDs,respAdj,raw,myK);

catRTs = cell(myK,8,2);
comboPairs = {[1,5],[2,6],[3,7],[4,8]};
for ik = 1:myK
    figure();
%     sp(2,2,1);
    myChans = find(sortIDs(:,myK)==ik);
    clustRTsPro = [];
    for ic = 1:length(myChans)
        for it = 1:8
            for ii = 1:2
                catRTs{ik,it,ii} = cat(1,catRTs{ik,it,ii},paRTs{myChans(ic),it,ii});
            end
        end
    end
    for i = 1:4
        yyaxis left
        sp(i) = subplot(2,2,i); hold on;
        [~,allBins] = ksdensity([catRTs{ik,comboPairs{i}(1),1};catRTs{ik,comboPairs{i}(1),2};catRTs{ik,comboPairs{i}(2),1};catRTs{ik,comboPairs{i}(2),2}],'bandwidth',20);
        [hcIn,hcBins] = ksdensity(catRTs{ik,comboPairs{i}(1),1},allBins,'bandwidth',20);
        heIn = ksdensity(catRTs{ik,comboPairs{i}(1),2},allBins,'bandwidth',20);
        hcOut = ksdensity(catRTs{ik,comboPairs{i}(2),1},allBins,'bandwidth',20);
        heOut = ksdensity(catRTs{ik,comboPairs{i}(2),2},allBins,'bandwidth',20);
        plot(allBins,hcIn,'color',[.2 .8 .2],'linewidth',3);
        plot(allBins,hcOut,'color',[.2 .8 .2],'linewidth',1);%'linestyle',':');
        plot(allBins,heIn,'color',[.8 .2 .2],'linewidth',3);
        plot(allBins,heOut,'color',[.8 .2 .2],'linewidth',1);%'linestyle',':');
        if ismember(i,[1,3])
            ylabel('RT PDF');
        end
        percErr = (sum(~isnan(catRTs{ik,comboPairs{i}(1),2}))+sum(~isnan(catRTs{ik,comboPairs{i}(2),2})))./(sum(~isnan(catRTs{ik,comboPairs{i}(1),2}))+sum(~isnan(catRTs{ik,comboPairs{i}(2),2}))+sum(~isnan(catRTs{ik,comboPairs{i}(1),1}))+sum(~isnan(catRTs{ik,comboPairs{i}(2),1})));
        yyaxis right
        plot([1000,1000],[0,percErr],'color',[.3 .3 .3],'linewidth',4); set(gca,'YLim',[0,1]);
        if ismember(i,[2,4])
            ylabel('Error Rate');
        end
    end
    set(sp,'XLim',[0,1000]);
    axes(sp(3)); xlabel('RT Stim-Aligned Trials');
    axes(sp(4)); xlabel('RT Sacc-Aligned Trials');
    suptitle(sprintf('Pro/Anti Cluster %d RTs',ik));
    clear sp
    
end

        
    
    
    
    
    
    
    
    

