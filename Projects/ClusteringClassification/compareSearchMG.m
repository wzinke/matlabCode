function compareSearchMG()

clear all; close all;

reload = 0;
reClust = 1;
areaCrit = {'FEF'};%,'F2','SEF','SC'};
recCrit = {'plex','tdt'};%{'rich','brad','plex','tdt'};
blWind = -300:-100;
searchWind = {[-50,200],[-150,100]};
doReward = 0;
minN = 10;
pAlph = .05;
saveLoad = 1;
doClip = 1;
matLoc = '~/Dropbox/Schall-Lab/dataMats';
saveDir = '~/Dropbox/Schall-Lab/Figures/fefMGtoSearchClusts';

%% Get a list of good variables
goodVars = whos; keepNames = {goodVars.name}; keepNames{end+1} = 'keepNames'; keepNames{end+1} = 'goodVars';

%% Start with MG
if reload
    [mgSDF,mgTimes,mgAreas,mgRecSys,mgRows,mgMonks,mgOtherTasks,mgSess,mgChans] = klPullAllSDFs('mg','-r',doReward,'rec',recCrit,'-a',areaCrit);
else
    % Load memory guided
    mgFiles = dir([matLoc,'/memGuide/*r',num2str(doReward),'.mat']);
    mgDates = cellfun(@(x) str2double(x(8:13)),{mgFiles.name});
    load([matLoc,'/memGuide/',mgFiles(mgDates==max(mgDates)).name]);
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
inSess = mgSess(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit));
inChans = mgChans(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit));

% Convert those SDFs to a cell that can be input to the clustering
inSDFs{1} = mgSDF{1}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
inSDFs{2} = mgSDF{2}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
if doReward
    inSDFs{3} = mgSDF{3}(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit),:);
end

% Get time cell
inTimes{1} = nanmean(cell2mat(mgTimes(:,1)),1);
inTimes{2} = nanmean(cell2mat(mgTimes(:,2)),1);
if doReward
    inTimes{3} = nanmean(cell2mat(mgTimes(:,3)),1);
end
mgParams = mgParams(ismember(mgAreas,areaCrit) & ismember(mgRecSys,recCrit));

% Make sure we're only using complete observations
catNorms = [];
for i = 1:length(inSDFs)
    catNorms = cat(2,catNorms,inSDFs{i}(:,ismember(inTimes{i},myEpocWinds{i})));
end
goodRows = sum(isfinite(catNorms),2) == size(catNorms,2);

inAreas = inAreas(goodRows);
inSess = inSess(goodRows);
inChans = inChans(goodRows);

for i = 1:length(inSDFs)
    inSDFs{i} = inSDFs{i}(goodRows,:);
end

% Re-assign and clear out
mgSDFs = inSDFs;
mgTimes = inTimes;
mgAreas = inAreas;
mgSess = inSess;
mgChans = inChans;
mgParams = mgParams(goodRows);

% Cluster
[mgIDs, ~, mgRaw, ~, mgLink] = klDoMetaClusteringv2(inSDFs,inTimes,'-e',myEpocs,'-ei',myEpocInds,'-er',myEpocWinds,'-r',0,'-mn',minN);

nNan = nan(1,size(mgIDs,2));
nBad = sum(isnan(mgIDs(:,1)));
nGood = sum(~isnan(mgIDs(:,1)));
for ik = 1:size(mgIDs,2)
    nNan(ik) = (sum(isnan(mgIDs(:,ik)))-nBad)./nGood;
end
mgK = find(nNan <= .1,1,'last');

% Now clear out
currVars = whos; currNames = {currVars.name}; currNames{end+1} = 'currNames'; currNames{end+1} = 'currVars'; currNames{end+1} = 'clearVars';
keepNames{end+1} = 'mgSDFs'; 
keepNames{end+1} = 'mgTimes'; 
keepNames{end+1} = 'mgAreas'; 
keepNames{end+1} = 'mgSess'; 
keepNames{end+1} = 'mgChans'; 
keepNames{end+1} = 'mgParams';
keepNames{end+1} = 'mgIDs';
keepNames{end+1} = 'mgRaw';
keepNames{end+1} = 'mgLink';
keepNames{end+1} = 'mgK';
keepNames{end+1} = 'clearVars';
clearVars = currNames(~ismember(currNames,keepNames));
cellfun(@clear,clearVars);

%% Prepare Search SDFs
if reload
    [searchSDF, searchTimes, searchAreas, searchRecSys, searchRow, searchMonks, searchOtherTasks, searchSess, searchChans] = klPullAllSDFs('search','-a',areaCrit,'-r',doReward,'-c',doClip);
else
    % Load search
    searchFiles = dir([matLoc,'/search/*r',num2str(doReward),'c',num2str(doClip),'.mat']);
    searchDates = cellfun(@(x) str2double(x(12:17)),{searchFiles.name});
    load([matLoc,'/search/',searchFiles(searchDates==max(searchDates)).name]);
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
myEpocWinds = {[-100:150],[-150:100],[-100:150],[-150:100]};
if doReward
    myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preRwd,earlyRwd,lateRwd,preVis,visTrans,visSust,preMov,postMov,nextVis,preRwd,earlyRwd,lateRwd};
    myEpocInds = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6];
    myEpocWinds = {[-200:300],[-300:200],[-100:100],[-200:300],[-300:200],[-100:100]};%myEpocWinds = cat(2,myEpocWinds,{[-100:100],[-100:100]});
end

inAreas = searchAreas(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit));
inRecSys = searchRecSys(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit));
inSess = searchSess(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit));
inChans = searchChans(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit));

% Convert those SDFs to a cell that can be input to the clustering
inSDFs{1} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),1),'UniformOutput',0));
inSDFs{2} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),2),'UniformOutput',0));
inSDFs{3} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),1),'UniformOutput',0));
inSDFs{4} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),2),'UniformOutput',0));
normSD = [1:4];
if doReward
    inSDFs{5} = inSDFs{4};
    inSDFs{4} = inSDFs{3};
    inSDFs{3} = cell2mat(cellfun(@(x) x(1,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),3),'UniformOutput',0));
    inSDFs{6} = cell2mat(cellfun(@(x) x(2,:),searchSDF(ismember(searchAreas,areaCrit) & ismember(searchRecSys,recCrit),3),'UniformOutput',0));
    normSD = cat(2,normSD,[5,6]);
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
end

inAreas = inAreas(goodRows);
inSess = inSess(goodRows);
inChans = inChans(goodRows);

% Re-assign and clear out
searchSDFs = inSDFs;
searchTimes = inTimes;
searchAreas = inAreas;
searchSess = inSess;
searchChans = inChans;

% Cluster
[searchIDs, ~, searchRaw, searchSumStruct, searchLink] = klDoMetaClusteringv2(searchSDFs,searchTimes,'-e',myEpocs,'-ei',myEpocInds,'-er',myEpocWinds,'-r',0,'sd',normSD,'-mn',minN);

nNan = nan(1,size(searchIDs,2));
nBad = sum(isnan(searchIDs(:,1)));
nGood = sum(~isnan(searchIDs(:,1)));
for ik = 1:size(searchIDs,2)
    nNan(ik) = (sum(isnan(searchIDs(:,ik)))-nBad)./nGood;
end
searchK = find(nNan <= .1,1,'last');

% Now clear out
currVars = whos; currNames = {currVars.name}; currNames{end+1} = 'currNames'; currNames{end+1} = 'currVars'; currNames{end+1} = 'clearVars';
keepNames{end+1} = 'searchSDFs'; 
keepNames{end+1} = 'searchTimes'; 
keepNames{end+1} = 'searchAreas'; 
keepNames{end+1} = 'searchSess'; 
keepNames{end+1} = 'searchChans';
keepNames{end+1} = 'searchIDs';
keepNames{end+1} = 'searchRaw';
keepNames{end+1} = 'searchLink';
keepNames{end+1} = 'searchK';
keepNames{end+1} = 'clearVars';
clearVars = currNames(~ismember(currNames,keepNames));
cellfun(@clear,clearVars);

%% Get units that are in both tasks - apply findMatch to the MG data
searchPair = cellfun(@(x,y) findMatch(x,y,searchSess,searchChans),mgSess,mgChans);
mgPair = cellfun(@(x,y) findMatch(x,y,mgSess,mgChans),searchSess,searchChans);

%% Normalize SDFs
normMG = klNormRespv2(mgSDFs,mgTimes,'ztr');
normSearch = klNormRespv2(searchSDFs,searchTimes,'ztr');

mgPairClusts = mgIDs(isfinite(searchPair),mgK);
uClusts = nunique(mgPairClusts);
svWind = searchWind{1}(1):searchWind{1}(2); smWind = searchWind{2}(1):searchWind{2}(2);
mgColors = jet(mgK);
for ic = 1:length(uClusts)
    figure(ic);
    sp(ic,1) = subplot(1,2,1); hold on;
    plot(searchTimes{1}(ismember(searchTimes{1},svWind)),normSearch{1}(searchPair(mgIDs(:,mgK)==uClusts(ic) & isfinite(searchPair)),ismember(searchTimes{1},svWind)),'color',[.8 .2 .2]);
    plot(searchTimes{3}(ismember(searchTimes{3},svWind)),normSearch{3}(searchPair(mgIDs(:,mgK)==uClusts(ic) & isfinite(searchPair)),ismember(searchTimes{3},svWind)),'color',[.2 .2 .8]);
    sp(ic,2) = subplot(1,2,2); hold on;
    plot(searchTimes{2}(ismember(searchTimes{2},smWind)),normSearch{2}(searchPair(mgIDs(:,mgK)==uClusts(ic) & isfinite(searchPair)),ismember(searchTimes{2},smWind)),'color',[.8 .2 .2]);
    plot(searchTimes{4}(ismember(searchTimes{4},smWind)),normSearch{4}(searchPair(mgIDs(:,mgK)==uClusts(ic) & isfinite(searchPair)),ismember(searchTimes{4},smWind)),'color',[.2 .2 .8]);
    suptitle(sprintf('MG Cluster %d - Individual Search Responses',uClusts(ic)));
    
    figure(100+ic);
    ap(ic,1) = subplot(1,2,1); hold on;
    plot(searchTimes{1}(ismember(searchTimes{1},svWind)),nanmean(normSearch{1}(searchPair(mgIDs(:,mgK)==uClusts(ic) & isfinite(searchPair)),ismember(searchTimes{1},svWind)),1),'color',[.8 .2 .2]);
    plot(searchTimes{3}(ismember(searchTimes{3},svWind)),nanmean(normSearch{3}(searchPair(mgIDs(:,mgK)==uClusts(ic) & isfinite(searchPair)),ismember(searchTimes{1},svWind)),1),'color',[.2 .2 .8]);
    ap(ic,2) = subplot(1,2,2); hold on;
    plot(searchTimes{2}(ismember(searchTimes{2},smWind)),nanmean(normSearch{2}(searchPair(mgIDs(:,mgK)==uClusts(ic) & isfinite(searchPair)),ismember(searchTimes{2},smWind)),1),'color',[.8 .2 .2]);
    plot(searchTimes{4}(ismember(searchTimes{4},smWind)),nanmean(normSearch{4}(searchPair(mgIDs(:,mgK)==uClusts(ic) & isfinite(searchPair)),ismember(searchTimes{4},smWind)),1),'color',[.2 .2 .8]);
    suptitle(sprintf('MG Cluster %d - Search Responses',uClusts(ic)));
    
    figure(200+ic);
    mp(ic,1) = subplot(1,2,1); hold on;
    plot(mgTimes{1},nanmean(normMG{1}(mgIDs(:,mgK)==uClusts(ic) & ~isnan(searchPair),:),1),'color',mgColors(ic,:));
    mp(ic,2) = subplot(1,2,2); hold on;
    plot(mgTimes{2},nanmean(normMG{2}(mgIDs(:,mgK)==uClusts(ic) & ~isnan(searchPair),:),1),'color',mgColors(ic,:));
    suptitle(sprintf('MG Cluster %d - MG Responses',uClusts(ic)));
    
    figure(300+ic);
    am(ic,1) = subplot(1,2,1); hold on;
    plot(mgTimes{1},nanmean(normMG{1}(mgIDs(:,mgK)==uClusts(ic),:),1),'color',mgColors(ic,:));
    am(ic,2) = subplot(1,2,2); hold on;
    plot(mgTimes{2},nanmean(normMG{2}(mgIDs(:,mgK)==uClusts(ic),:),1),'color',mgColors(ic,:));
    suptitle(sprintf('MG Cluster %d - All MG Responses',uClusts(ic)));
    
end
set(sp(:,1),'XLim',searchWind{1},'tickdir','out','box','off','ticklength',get(sp(1,1),'ticklength').*3);
set(sp(:,2),'XLim',searchWind{2},'yaxisloc','right','tickdir','out','box','off','ticklength',get(sp(1,2),'ticklength').*3);
set(ap(:,1),'XLim',searchWind{1},'tickdir','out','box','off','ticklength',get(ap(1,1),'ticklength').*3);
set(ap(:,2),'XLim',searchWind{2},'yaxisloc','right','tickdir','out','box','off','ticklength',get(ap(1,2),'ticklength').*3);
set(mp(:,1),'XLim',[-200,300],'tickdir','out','box','off','ticklength',get(mp(1,1),'ticklength').*3);
set(mp(:,2),'XLim',[-300,200],'yaxisloc','right','tickdir','out','box','off','ticklength',get(mp(1,2),'ticklength').*3);
set(am(:,1),'XLim',[-200,300],'tickdir','out','box','off','ticklength',get(am(1,1),'ticklength').*3);
set(am(:,2),'XLim',[-300,200],'yaxisloc','right','tickdir','out','box','off','ticklength',get(am(1,2),'ticklength').*3);
linkaxes(sp,'y');
linkaxes(ap,'y');
linkaxes(mp,'y');
linkaxes(am,'y');

strs = {'01','02','03','04','05','06','07','08','09','10'};
for ic = 1:length(uClusts)
    saveas(ic,[saveDir,'/mgClusts/searchResps-mgClust-indivs',strs{ic},'.fig']);
    saveas(ic+100,[saveDir,'/mgClusts/searchResps-mgClust',strs{ic},'.fig']);
    saveas(ic+200,[saveDir,'/mgClusts/mgResps-mgClust',strs{ic},'.fig']);
    saveas(ic+300,[saveDir,'/mgClusts/mgRespsAll-mgClust',strs{ic},'.fig']);
end
close all;

%% Run some stats 

% Check Target Selection in Search
allTST = getDiff(normSearch{1},normSearch{3},searchTimes{1});
% [p,t,stats] = anovan(allTST(searchPair(isfinite(searchPair))),mgIDs(isfinite(searchPair),mgK));

% Extract mg spiking parameters
mgTuneWidth = reshape([mgParams.sigma],2,length(mgParams))';
mgTuneDir = reshape([mgParams.theta],2,length(mgParams))';
mgTuneAmp = reshape([mgParams.alpha],2,length(mgParams))';
mgBaseline = [mgParams.blRate]';
tmpV = [mgParams.vGOF]';
mgVFitSSE = [tmpV.sse]';
mgVFitR2 = [tmpV.rsquare]';
tmpM = [mgParams.mGOF]';
mgMFitSSE = [tmpM.sse]';
mgMFitR2 = [tmpM.rsquare]';
mgCV = [mgParams.CV]';
mgCV2 = [mgParams.CV2]';
mgLV = [mgParams.LV]';
mgLVR = [mgParams.LVR]';
spkWidths = [mgParams.spkWidth]';

anovaVars = {'mgTuneWidth(mgVFitR2 >= .5,1)','mgTuneWidth(mgMFitR2 >= .5,2)','mgTuneDir(mgVFitR2 >= .5,1)','mgTuneDir(mgMFitR2 >= .5,2)','mgTuneAmp(mgVFitR2 >= .5,1)','mgTuneAmp(mgMFitR2 >= .5,2)',...
    'mgBaseline','mgCV','mgCV2','mgLV','mgLVR','allTST(searchPair(isfinite(searchPair)))','spkWidths'};
anovaSubset = {'mgVFitR2 >= .5','mgMFitR2 >= .5','mgVFitR2 >= .5','mgMFitR2 >= .5','mgVFitR2 >= .5','mgMFitR2 >= .5',...
    ':',':',':',':',':','isfinite(searchPair)',':'};
anovaYStr = {'Visual Tuning Width (deg.)','Movement Tuning Width (deg.)','Visual RF (deg.)','Movement Field (deg).','Visual Max (sp/s)','Movement Max (sp/s)',...
    'Baseline FR (sp/s)','CV','CV2','LV','LVR','Target Selection Time (ms)','Spike Width'};
anovaSaveStr = {'vTuneWd','mTuneWd','vTuneDir','mTuneDir','vRF','mRF','vMax','mMax','blRate','CV','CV2','LV','LVR','TST','spkWidth'};

for iv = 1:length(anovaVars)
%     eval(['[p(iv),t{iv},stats(iv)] = anovan(',anovaVars{iv},',mgIDs(',anovaSubset{iv},',mgK),''display'',''off'');']);
    eval(['[p(iv),t{iv},stats(iv)] = kruskalwallis(',anovaVars{iv},',mgIDs(',anovaSubset{iv},',mgK),''off'');']);    
end

sigComps = find(p < pAlph);
for ip = 1:length(sigComps)
    % Make the scatter plot
    
    figure(); %hold on;
    eval(['scatter(mgIDs(',anovaSubset{sigComps(ip)},',mgK),',anovaVars{sigComps(ip)},',[],''k'');']);
    hold on;
    % Find the post-hoc comparisons
    mcMat = multcompare(stats(sigComps(ip)),'display','off');
    postHocSig = find(mcMat(:,end) <= pAlph);
    if ~isempty(postHocSig)
        % Get default Y axis values
        currY = get(gca,'YLim');
        % Let's let the post-hoc indicators span 10% at the top
        for ii = 1:length(postHocSig)
            plot([mcMat(postHocSig(ii),1),mcMat(postHocSig(ii),2)],ones(1,2).*currY(2)+((.1.*ii*range(currY))/length(postHocSig)),'color','k');
        end
        set(gca,'YLim',[currY(1),currY(2)+.12*range(currY)]);
    end
    set(gca,'tickdir','out','ticklength',get(gca,'ticklength').*3,'box','off','XLim',[0,mgK+1]);
    xlabel('MG Cluster');
    ylabel(anovaYStr{sigComps(ip)});
    saveas(gcf,[saveDir,'/stats/aov-',anovaSaveStr{sigComps(ip)},'.fig']);
    clear mcMat
end

clear p

% Now let's do the comparison of clusters
close all; 
valMG = cellfun(@(x) x(isfinite(searchPair),:),normMG,'UniformOutput',0);
valSearch = cellfun(@(x) x(searchPair(isfinite(searchPair)),:),normSearch,'UniformOutput',0);
valMGClusts = mgIDs(isfinite(searchPair),mgK);
valSearchClusts = searchIDs(searchPair(isfinite(searchPair)),searchK);

[p,x2,~,exp,in,~,congingY,contingX] = klConting(valMGClusts,valSearchClusts);
rDiff = (in-exp)./exp;
for i = 1:1000
    [~,~,~,expRand,inRand] = klConting(klShuffle(valMGClusts),klShuffle(valSearchClusts));
    randDiff(:,:,i) = (inRand-expRand)./expRand;
end
zDiff = (rDiff-nanmean(randDiff,3))./nanstd(randDiff,[],3);
figure(1000);
imagesc(zDiff);
cb = colorbar;
caxis([-4,4]);
colormap('jet');
set(gca,'box','off','XTick',[],'YTick',[]);
set(cb,'tickdir','out','ticklength',get(cb,'ticklength').*3,'Ticks',[],'TickLabels',[]);
saveas(gcf,[saveDir,'/mgSearchComparison-Heatmap.fig']);

% Let's do the same as a percentage of a given MG cluster belonging to a
% search cluster
percMatM = in./repmat(sum(in,2),1,size(in,2));
figure(1001);
imagesc(percMatM);
cb2 = colorbar;
caxis([0,1]);
colormap('hot');
set(gca,'box','off','XTick',[],'YTick',[]);
set(cb2,'tickdir','out','ticklength',get(cb2,'ticklength').*3,'Ticks',[],'TickLabels',[]);
saveas(gcf,[saveDir,'/mgSearchComparison-PercMGClust-Heatmap.fig']);

% Let's do the same as a percentage of a given search cluster belonging to a
% MG cluster
percMatS = in./repmat(sum(in,1),size(in,1),1);
figure(1002);
imagesc(percMatS);
cb3 = colorbar;
caxis([0,1]);
colormap('hot');
set(gca,'box','off','XTick',[],'YTick',[]);
set(cb3,'tickdir','out','ticklength',get(cb3,'ticklength').*3,'Ticks',[],'TickLabels',[]);
saveas(gcf,[saveDir,'/mgSearchComparison-PercSearchClust-Heatmap.fig']);




% Now make search clusters
uClusts = nunique(valSearchClusts);
for ic = 1:length(uClusts)
%     figure(ic); hold on;
%     figure(ic+100); hold on;
    for ip = 1:2
        figure(ic);
        smp(ic,ip) = subplot(1,2,ip); hold on;
        plot(searchTimes{ip},nanmean(valSearch{ip}(valSearchClusts==uClusts(ic),:),1),'color','r');
        plot(searchTimes{ip+2},nanmean(valSearch{ip+2}(valSearchClusts==uClusts(ic),:),1),'color','b');
        
        figure(ic+100);
        sap(ic,ip) = subplot(1,2,ip); hold on;
        plot(searchTimes{ip},valSearch{ip}(valSearchClusts==uClusts(ic),:),'color','r');
        plot(searchTimes{ip+2},valSearch{ip+2}(valSearchClusts==uClusts(ic),:),'color','b');
        
        figure(ic+200);
        saa(ic,ip) = subplot(1,2,ip); hold on;
        plot(searchTimes{ip},normSearch{ip}(searchIDs(:,searchK)==uClusts(ic),:),'color','r');
        plot(searchTimes{ip+2},normSearch{ip+2}(searchIDs(:,searchK)==uClusts(ic),:),'color','b');
        
        figure(ic+300);
        sma(ic,ip) = subplot(1,2,ip); hold on;
        plot(searchTimes{ip},nanmean(normSearch{ip}(searchIDs(:,searchK)==uClusts(ic),:),1),'color','r');
        plot(searchTimes{ip+2},nanmean(normSearch{ip+2}(searchIDs(:,searchK)==uClusts(ic),:),1),'color','b');
    end
    figure(ic);
    suptitle(['Search Responses - Search Cluster ',num2str(ic)]);
    figure(ic+100);
    suptitle(['Indiv Search Responses - Search Cluster ',num2str(ic)]);
    figure(ic+200);
    suptitle(['All Indiv Search Responses - Search Cluster ',num2str(ic)]);
    figure(ic+300);
    suptitle(['All Search Responses - Search Cluster ',num2str(ic)]);
    
end
set(smp(:,1),'XLim',searchWind{1},'tickdir','out','box','off','ticklength',get(smp(1,1),'ticklength').*3);
set(smp(:,2),'XLim',searchWind{2},'yaxisloc','right','tickdir','out','box','off','ticklength',get(smp(1,2),'ticklength').*3);
set(sap(:,1),'XLim',searchWind{1},'tickdir','out','box','off','ticklength',get(sap(1,1),'ticklength').*3);
set(sap(:,2),'XLim',searchWind{2},'yaxisloc','right','tickdir','out','box','off','ticklength',get(sap(1,2),'ticklength').*3);
set(saa(:,1),'XLim',searchWind{1},'tickdir','out','box','off','ticklength',get(saa(1,1),'ticklength').*3);
set(saa(:,2),'XLim',searchWind{2},'yaxisloc','right','tickdir','out','box','off','ticklength',get(saa(1,2),'ticklength').*3);
set(sma(:,1),'XLim',searchWind{1},'tickdir','out','box','off','ticklength',get(sma(1,1),'ticklength').*3);
set(sma(:,2),'XLim',searchWind{2},'yaxisloc','right','tickdir','out','box','off','ticklength',get(sma(1,2),'ticklength').*3);
linkaxes(smp,'y');
linkaxes(sap,'y');
linkaxes(saa,'y');
linkaxes(sma,'y');

for ic = 1:length(uClusts)
    saveas(ic,[saveDir,'/searchClusts/searchResps-searchClust-indivs',strs{ic},'.fig']);
    saveas(ic+100,[saveDir,'/searchClusts/searchResps-searchClust',strs{ic},'.fig']);
    saveas(ic+200,[saveDir,'/searchClusts/searchResps-searchClust-allIndivs',strs{ic},'.fig']);
    saveas(ic+300,[saveDir,'/searchClusts/searchResps-searchClust-all',strs{ic},'.fig']);
end


keyboard

end

function ind = findMatch(sourceSess,sourceChan,matchSess,matchChan)

ind = find(ismember(matchSess,sourceSess) & ismember(matchChan,sourceChan));
if isempty(ind), ind = nan; end

end

function [diffT, diffSDF] = getDiff(val1,val2,times)
    % Get difference and z-score
    diffSDF = val1 - val2;
    blMean = nanmean(diffSDF(:,times >= -200 & times <= -100),2);
    blStd = nanstd(diffSDF(:,times >= -200 & times <= -100),[],2);
    zSDF = (diffSDF-repmat(blMean,1,size(diffSDF,2)))./repmat(blStd,1,size(diffSDF,2));
    
    % Get the first time that zSDF >= 2 for more than 15ms
    overFive = klGetConsecutive(zSDF >= 6);
    overTwo = klGetConsecutive(zSDF >= 2);
    [overR,overC] = find(overFive >= 10 & times >= 0);
    diffT = nan(size(val1,1),1);
    for ii = 1:length(diffT)
        if any(overR==ii)
            diffT(ii) = times(min(overC(overR==ii & overC >= find(times == 0))));
            if diffT(ii) > 200
                diffT(ii) = nan;
            end
        end
    end
     
end
