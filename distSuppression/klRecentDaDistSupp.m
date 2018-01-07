function klRecentDaDistSupp(varargin)
clearvars; close all;

tebaDir = [tebaMount,'Users/Kaleb/dataProcessed'];
areaCrit = {'FEF'};
blWind = -300:-100;
tstWind = 100:150;
rePull = 0;
cutChecks = 0;

% Manually select Da rows
daGoodRows = [80 144 149 152 168 169 205 206 208 217 222 223 225 240 241 242 296 297 299 305 309 312 319 321 323 367 372 436 533 535 537 540 541 542 554 585 586 602 657 669 672 673 674 682 683 685 686]';
checkRFs = [80 168 217 240 299 321 323 537 540 541 602 657];
cutAnyway = [208 223 309 312 319 533 535 542 554 672 683 80 152 168 602];
daGoodRows(ismember(daGoodRows,cutAnyway)) = [];
if cutChecks
    daGoodRows(ismember(daGoodRows,checkRFs)) = [];
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
myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis,preVis,visTrans,visSust,preMov,postMov,nextVis,preVis,visTrans,visSust,preMov,postMov,nextVis};
myEpocInds = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6];
myEpocWinds = {[-200:300],[-300:200],[-200:300],[-300:200],[-200:300],[-300:200]};

% Get all capture SDFs
if rePull
    [capSDF,capSDFTimes,capAreas,capRows,capMonks,capVis,capTST,capDST,capPerf,capRF,capRT,capTSTA] = klPullKilosortSDFs('cap','-a',{'FEF'},'-k','psp','-w',20,'rows',daGoodRows);
    capAreas = cell(length(capDST),1);
    [capAreas{1:length(capAreas)}] = deal('FEF');
    save('./capSDFs-allNewP.mat','cap*');
else
    load('capSDFs-allNewP.mat');
end

capSizes = cellfun(@size,capSDF,'UniformOutput',0);
hasEmpty = ~logical(sum(cellfun(@(x) any(x==0),capSizes(:,1:2)),2));

capVars = whos('cap*');
for ic = 1:length(capVars)
    eval([capVars(ic).name,'=',capVars(ic).name,'(hasEmpty,:);']);
end

% Convert those SDFs to a cell that can be input to the clustering
inSDFs{1} = cell2mat(cellfun(@(x) x(3,:),capSDF(ismember(capAreas,areaCrit),1),'UniformOutput',0));
inSDFs{2} = cell2mat(cellfun(@(x) x(3,:),capSDF(ismember(capAreas,areaCrit),2),'UniformOutput',0));
inSDFs{3} = cell2mat(cellfun(@(x) x(5,:),capSDF(ismember(capAreas,areaCrit),1),'UniformOutput',0));
inSDFs{4} = cell2mat(cellfun(@(x) x(5,:),capSDF(ismember(capAreas,areaCrit),2),'UniformOutput',0));
inSDFs{5} = cell2mat(cellfun(@(x) x(4,:),capSDF(ismember(capAreas,areaCrit),1),'UniformOutput',0));
inSDFs{6} = cell2mat(cellfun(@(x) x(4,:),capSDF(ismember(capAreas,areaCrit),2),'UniformOutput',0));
normSD = [1:6];

plotSDFs{1} = cell2mat(cellfun(@(x) x(3,:),capSDF(ismember(capAreas,areaCrit),1),'UniformOutput',0));
plotSDFs{2} = cell2mat(cellfun(@(x) x(3,:),capSDF(ismember(capAreas,areaCrit),2),'UniformOutput',0));
plotSDFs{3} = cell2mat(cellfun(@(x) x(5,:),capSDF(ismember(capAreas,areaCrit),1),'UniformOutput',0));
plotSDFs{4} = cell2mat(cellfun(@(x) x(5,:),capSDF(ismember(capAreas,areaCrit),2),'UniformOutput',0));
plotSDFs{5} = cell2mat(cellfun(@(x) x(4,:),capSDF(ismember(capAreas,areaCrit),1),'UniformOutput',0));
plotSDFs{6} = cell2mat(cellfun(@(x) x(4,:),capSDF(ismember(capAreas,areaCrit),2),'UniformOutput',0));

% Get time cell
inTimes{1} = nanmean(cell2mat(capSDFTimes(:,1)),1);
inTimes{2} = nanmean(cell2mat(capSDFTimes(:,2)),1);
inTimes{3} = inTimes{1};
inTimes{4} = inTimes{2};
inTimes{5} = inTimes{1};
inTimes{6} = inTimes{2};

catNorms = [];
for i = 1:length(inSDFs)
    catNorms = cat(2,catNorms,inSDFs{i}(:,ismember(inTimes{i},myEpocWinds{i})));
end
goodRows = sum(isfinite(catNorms),2) == size(catNorms,2);

for i = 1:length(inSDFs)
    inSDFs{i} = inSDFs{i}(goodRows,:);
    plotSDFs{i} = plotSDFs{i}(goodRows,:);
end
% capRows = capRows(goodRows);
for ic = 1:length(capVars)
    eval([capVars(ic).name,'=',capVars(ic).name,'(goodRows,:);']);
end


% Use capRows to get the appropriate sessions
[num,text,all]=xlsread('kiloSortBookKeeping.xlsx');
uSessProbe = unique(all(capRows,1));
uSess = cellfun(@(x) x(2:(strfind(x,'_probe')-1)),uSessProbe,'UniformOutput',0);
uSess = unique(uSess);

whichSess = nan(length(capRows),1);
unitAcc = nan(length(capRows),1);

for is = 1:length(uSess)
    clear Task
    % Identify sessions to units
    whichSess(ismember(all(capRows,1),uSessProbe{is})) = is;
    
    % Load in Behav.mat
    load([tebaDir,filesep,uSess{is},'/Behav.mat']);
    
    % Get Percent correct for capture vs no capture
    percCap = sum(strcmpi(Task.TaskType,'Cap') & Task.Correct==1 & Task.Singleton==2222)./sum(strcmpi(Task.TaskType,'Cap') & (Task.Error == 5 | Task.Correct==1) & Task.Singleton==2222);
    percTrad = sum(strcmpi(Task.TaskType,'Cap') & Task.Correct==1 & Task.Singleton==1111)./sum(strcmpi(Task.TaskType,'Cap') & (Task.Error == 5 | Task.Correct==1) & Task.Singleton==1111);
    
    capRTs = nanmean(Task.SRT(strcmpi(Task.TaskType,'Cap') & Task.Correct==1 & Task.Singleton==2222));
    tradRT = nanmean(Task.SRT(strcmpi(Task.TaskType,'Cap') & Task.Correct==1 & Task.Singleton==1111));
    
    percCorrect(is,:) = [percCap,percTrad];
    sessRTs(is,:) = [capRTs,tradRT];
    unitAcc(ismember(all(capRows,1),uSessProbe{is})) = percCorrect(is,1);
    
end

tInMinusOutV = inSDFs{1}-inSDFs{5};
tInMinusOutVZ = (tInMinusOutV-repmat(nanmean(tInMinusOutV(:,ismember(inTimes{1},blWind)),2),1,size(tInMinusOutV,2)))./repmat(nanstd(tInMinusOutV(:,ismember(inTimes{1},blWind)),[],2),1,size(tInMinusOutV,2));

% hasTST = any(klGetConsecutive(abs(tInMinusOutVZ) >= 5) >= 20,2);
sigSel=klGetConsecutive(abs(tInMinusOutVZ) >= 5);
allTST = nan(size(tInMinusOutVZ,1),1);
for i = 1:size(sigSel)
    if any(sigSel(i,:) >= 20 & inTimes{1} >= 0)
        allTST(i) = inTimes{1}(find(sigSel(i,:) >= 20 & inTimes{1} >= 0,1));
    end
end

% goodUnits=allTST < 200 & capVis;% & unitAcc >= .85;
% goodUnits = ~isnan(capTST) | ~isnan(capDST);
normResp = klNormRespv2(inSDFs,inTimes,'ztrbl','-r',myEpocWinds);
tsts = getDiff(normResp{1},normResp{5},inTimes{1});
dsts = getDiff(normResp{5},normResp{3},inTimes{1});

% figure(); hold on;
% plot(inTimes{1},nanmean(normResp{1}(goodUnits,:),1),'color',[.8 .2 .2]);
% plot(inTimes{3},nanmean(normResp{3}(goodUnits,:),1),'color',[.2 .8 .2]);
% plot(inTimes{5},nanmean(normResp{5}(goodUnits,:),1),'color',[.2 .2 .8]);
% title('Da-TDT','fontsize',18);
% xlabel('Time From Stimulus (ms)','fontsize',14);
% ylabel('Normalized Firing Rate (Z)','fontsize',14);
% set(gca,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',12,'XLim',[-100,300]);
% 
tDiff = nanmean(normResp{1}(:,ismember(inTimes{1},tstWind)),2)-nanmean(normResp{3}(:,ismember(inTimes{1},tstWind)),2);
dDiff = nanmean(normResp{5}(:,ismember(inTimes{1},tstWind)),2)-nanmean(normResp{3}(:,ismember(inTimes{1},tstWind)),2);

perfDiff = capPerf(:,2)-capPerf(:,1);
rtDiff = capRT(:,2) - capRT(:,1);
goodUnits = (~isnan(tsts) | ~isnan(dsts));

figure();
scatter(perfDiff(goodUnits),dDiff(goodUnits)); hold on;
[r,p] = corr(perfDiff(goodUnits),dDiff(goodUnits));
b=regress(dDiff(goodUnits),[perfDiff(goodUnits),ones(length(perfDiff(goodUnits)),1)]);
xLims = get(gca,'XLim');
plot(xLims,(xLims.*b(1))+b(2),'color','k');
ylabel('Nonsalient SDF - Salient SDF');
xlabel('Nonsalient Correct - Salient Correct');
title(sprintf('DaNew-r=%.3f-p=%.3f',r,p));
saveas(gcf,'./Figs/Da-New/DaNew-PercCorrectCorr.png');
% 
figure(); hold on;
% pltMeanStd(inTimes{1},nanmean(normResp{1}(goodUnits,:),1),nanstd(normResp{1}(goodUnits,:),[],1)./sqrt(sum(goodUnits)),'color',[.8 .2 .2]);
% pltMeanStd(inTimes{3},nanmean(normResp{3}(goodUnits,:),1),nanstd(normResp{3}(goodUnits,:),[],1)./sqrt(sum(goodUnits)),'color',[.2 .8 .2]);
% pltMeanStd(inTimes{5},nanmean(normResp{5}(goodUnits,:),1),nanstd(normResp{5}(goodUnits,:),[],1)./sqrt(sum(goodUnits)),'color',[.2 .2 .8]);
plot(inTimes{1},nanmean(normResp{1}(goodUnits,:),1),'color','k','linewidth',3,'linestyle','-');
plot(inTimes{3},nanmean(normResp{3}(goodUnits,:),1),'color','r','linewidth',1,'linestyle','-');
plot(inTimes{5},nanmean(normResp{5}(goodUnits,:),1),'color','k','linewidth',1,'linestyle','-');
plot([-50 -50],[1-nanmean(nanstd(normResp{1}(goodUnits,:),[],1)./sqrt(size(normResp{1}(goodUnits,:),1)),2),1+nanmean(nanstd(normResp{1}(goodUnits,:),[],1)./sqrt(size(normResp{1}(goodUnits,:),1)),2)],'color','k','linewidth',3);
plot([-25 -25],[1-nanmean(nanstd(normResp{3}(goodUnits,:),[],1)./sqrt(size(normResp{3}(goodUnits,:),1)),2),1+nanmean(nanstd(normResp{3}(goodUnits,:),[],1)./sqrt(size(normResp{3}(goodUnits,:),1)),2)],'color','r','linewidth',1);
plot([0 0],[1-nanmean(nanstd(normResp{5}(goodUnits,:),[],1)./sqrt(size(normResp{5}(goodUnits,:),1)),2),1+nanmean(nanstd(normResp{5}(goodUnits,:),[],1)./sqrt(size(normResp{5}(goodUnits,:),1)),2)],'color','k','linewidth',1);
title('Da-TDT','fontsize',18);
xlabel('Time From Stimulus (ms)','fontsize',14);
ylabel('Normalized Firing Rate (Z)','fontsize',14);
set(gca,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',12,'XLim',[-100,300]);
saveas(gcf,'./Figs/Da-New/DaNew-MnSDF.fig');

figure();
scatter(rtDiff,dDiff); hold on;
[r,p] = corr(rtDiff,dDiff);
b=regress(dDiff,[rtDiff,ones(length(rtDiff),1)]);
xLims = get(gca,'XLim');
plot(xLims,(xLims.*b(1))+b(2),'color','k');
ylabel('Nonsalient SDF - Salient SDF');
xlabel('Nonsalient RT - Salient RT');
title(sprintf('All-r=%.3f-p=%.3f',r,p));
saveas(gcf,'./Figs/Da-New/DaNew-RTCorr.png');
% saveas(gcf,['./Figs/All-',normType,'-RTCorr.png']);
% for im = 1:length(uMonks)
%     figure(); hold on;
%     scatter(rtDiff(ismember(goodMonks(goodUnits),uMonks{im})),dDiff(ismember(goodMonks(goodUnits),uMonks{im})));
%     [r,p] = corr(rtDiff(ismember(goodMonks(goodUnits),uMonks{im})),dDiff(ismember(goodMonks(goodUnits),uMonks{im})));
%     b=regress(dDiff(ismember(goodMonks(goodUnits),uMonks{im})),[rtDiff(ismember(goodMonks(goodUnits),uMonks{im})),ones(length(rtDiff(ismember(goodMonks(goodUnits),uMonks{im}))),1)]);
%     xLims = get(gca,'XLim');
%     plot(xLims,(xLims.*b(1))+b(2),'color','k');
%     ylabel('Nonsalient SDF - Salient SDF');
%     xlabel('Nonsalient RT - Salient RT');
%     title(sprintf('%s-r=%.3f-p=%.3f',uMonks{im},r,p));
% %     saveas(gcf,['./Figs/',uMonks{im},'-',normType,'-RTCorr.png']);
% end
keyboard

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



    
    
