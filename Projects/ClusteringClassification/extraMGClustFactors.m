function extraMGClustFactors()

clearvars; close all;

reload = 0;
reClust = 1;
areaCrit = {'FEF'};%,'F2','SEF','SC'};
recCrit = {'plex','tdt'};%{'rich','brad','plex','tdt'};
blWind = -300:-100;
searchWind = {[-50,200],[-150,100]};
doReward = 0;
minN = 10;
pAlph = .05;
pTrend = .1;
saveLoad = 1;
doClip = 1;
matLoc = [mlRoot,'Dropbox/Schall-Lab/dataMats'];
saveDir = [mlRoot,'Dropbox/Schall-Lab/Figures/fefMGtoSearchClusts'];

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

% Reorder according to paper clusters
reOrd = [5,2,4,10,8,9,1,3,6,7];
metaOrder = nan(size(mgIDs));
for i = 1:10
    metaOrder(mgIDs(:,mgK)==reOrd(i),mgK) = i;
end

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
mgFano = [mgParams.Fano]';
spkWidths = [mgParams.spkWidth]';

% Cut out silly values for fitting
mgTuneWidth(abs(mgTuneWidth) >= 180) = nan;
mgTuneAmp(abs(mgTuneAmp) >= 1000) = nan;


anovaVars = {'mgTuneWidth','mgTuneWidth','mgTuneDir','mgTuneDir','mgTuneAmp','mgTuneAmp',...
    'mgBaseline','mgCV','mgCV2','mgLV','mgLVR','mgFano','spkWidths'};
varSubset = {'mgVFitR2 >= .5','mgMFitR2 >= .5','mgVFitR2 >= .5','mgMFitR2 >= .5','mgVFitR2 >= .5','mgMFitR2 >= .5',...
    'true(length(mgBaseline),1)','true(length(mgCV),1)','true(length(mgCV2),1)','true(length(mgLV),1)','true(length(mgLVR),1)','true(length(mgFano),1)','true(length(spkWidths),1)'};
varCol = {',1',',2',',1',',2',',1',',2',[],[],[],[],[],[],[],[]};
clustSubset = {'mgVFitR2 >= .5','mgMFitR2 >= .5','mgVFitR2 >= .5','mgMFitR2 >= .5','mgVFitR2 >= .5','mgMFitR2 >= .5',...
    'true(length(mgBaseline),1)','true(length(mgCV),1)','true(length(mgCV2),1)','true(length(mgLV),1)','true(length(mgLVR),1)','true(length(mgFano),1)','true(length(spkWidths),1)'};
anovaYStr = {'Visual Tuning Width (deg.)','Movement Tuning Width (deg.)','Visual RF (deg.)','Movement Field (deg).','Visual Max (sp/s)','Movement Max (sp/s)',...
    'Baseline FR (sp/s)','CV','CV2','LV','LVR','Fano Factor','Spike Width'};
anovaSaveStr = {'vTuneWd','mTuneWd','vRF','mRF','vMax','mMax','blRate','CV','CV2','LV','LVR','Fano','spkWidth'};

for iv = 1:length(anovaVars)
%     eval(['[p(iv),t{iv},stats(iv)] = anovan(',anovaVars{iv},',metaOrder(',anovaSubset{iv},',mgK),''display'',''off'');']);
    eval(['[p(iv),t{iv},stats(iv)] = kruskalwallis(',anovaVars{iv},'(',varSubset{iv},varCol{iv},')',',metaOrder(',clustSubset{iv},',mgK),''off'');']);    
end

sigComps = find(p < pTrend);
for ip = 1:length(sigComps)
    % Make the scatter plot
    
    figure(); %hold on;
    eval(['scatter(metaOrder(',clustSubset{sigComps(ip)},',mgK),',anovaVars{sigComps(ip)},'(',varSubset{sigComps(ip)},varCol{sigComps(ip)},')',',[],''k'');']);
    hold on;
    % Find the post-hoc comparisons
    mcMat = multcompare(stats(sigComps(ip)),'display','off');
    postHocSig = find(mcMat(:,end) <= pTrend);
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
    title(sprintf('Omnibus %s - p=%.3f',anovaSaveStr{sigComps(ip)},p(sigComps(ip))));
    saveas(gcf,[saveDir,'/stats/aov-',anovaSaveStr{sigComps(ip)},'.fig']);
    clear mcMat postHocSig
end

clear p t stats sigComps
% close all;

% Now do some specific comparisons
% vmClusts = [1,4,9,10];
% mClusts = [3,6];
% vClusts = [2,5];
% clustsToCheck = {[1,4,9,10],[3,6],[2,5]};
clustsToCheck = {[3,4,6,7],[8,9],[1,2]};
clustAbbr = {'VM','M','V'};
for is = 1:length(clustsToCheck)
    clear sigComps mcMat postHocSig
    theseClusts = clustsToCheck{is};
    for iv = 1:length(anovaVars)
%         eval(['[p(is,iv),t{is,iv},stats.(clustAbbr{is})(iv)] = anovan(',anovaVars{iv},'(',varSubset{iv},'& ismember(metaOrder(:,mgK),theseClusts)',varCol{iv},')',',metaOrder(',clustSubset{iv},'& ismember(metaOrder(:,mgK),theseClusts),mgK),''display'',''off'');']);
        eval(['[p(is,iv),t{is,iv},stats.(clustAbbr{is})(iv)] = kruskalwallis(',anovaVars{iv},'(',varSubset{iv},'& ismember(metaOrder(:,mgK),theseClusts)',varCol{iv},')',',metaOrder(',clustSubset{iv},'& ismember(metaOrder(:,mgK),theseClusts),mgK),''off'');']);
    end
    
    sigComps = find(p(is,:) < pTrend);
    for ii = 1:length(sigComps)
        figure();
        eval(['scatter(metaOrder(ismember(metaOrder(:,mgK),theseClusts),mgK),',anovaVars{sigComps(ii)},'(ismember(metaOrder(:,mgK),theseClusts)',varCol{sigComps(ii)},'),[],''k'');']);
        hold on;
        mcMat = multcompare(stats.(clustAbbr{is})(sigComps(ii)),'display','off');
        postHocSig = find(mcMat(:,end) <= pTrend);
        currY = get(gca,'YLim');
        for iii = 1:length(postHocSig)
            plot([theseClusts(mcMat(postHocSig(iii),1)),theseClusts(mcMat(postHocSig(iii),2))],ones(1,2).*currY(2)+((.1.*iii*range(currY))/length(postHocSig)),'color','k');
        end
        set(gca,'YLim',[currY(1),currY(2)+.12*range(currY)]);
        set(gca,'tickdir','out','ticklength',get(gca,'ticklength').*3,'box','off','XLim',[min(theseClusts)-.5,max(theseClusts)+.5]);
        xlabel('MG Cluster');
        ylabel(anovaYStr{sigComps(ii)});
        title(sprintf('%s %s - p=%.3f',clustAbbr{is},anovaSaveStr{sigComps(ii)},p(is,sigComps(ii))));
        saveas(gcf,[saveDir,'/stats/aov-',clustAbbr{is},'-',anovaSaveStr{sigComps(ii)},'.fig']);
        clear mcMat postHocSig
            
    end
    
    
end