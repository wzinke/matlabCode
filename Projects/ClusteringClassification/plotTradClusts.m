
%% Get visual/movement identities (B&G-ish method)
checkVis = @(vSDF,vTimes) any(klGetConsecutive(abs(vSDF(:,vTimes >= 0 & vTimes <= 150)) > 6) >= 2);
checkMov = @(mSDF,mTimes) any(klGetConsecutive(abs(mSDF(:,mTimes <= 0 & mTimes >= -150)) > 6) >= 2);
[normVis, normMov] = klNormResp(goodSDF{1},allTimeCell{1},goodSDF{2},allTimeCell{2},'zbl');
for ii = 1:size(normVis,1),
isVis = checkVis(normVis(ii,:),allTimeCell{1});
[movCorr(ii),movP(ii)] = corr(normMov(ii,allTimeCell{2} >= -20 & allTimeCell{2} <= 0)',allTimeCell{2}(allTimeCell{2} >= -20 & allTimeCell{2} <= 0)');
isMov = checkMov(normMov(ii,:),allTimeCell{2}) && (movP(ii) < .05) && (movCorr(ii) > 0);
if isVis && isMov,
myTypes{ii} = 'vismov';
elseif isVis && ~isMov,
myTypes{ii} = 'vis';
elseif ~isVis && isMov,
myTypes{ii} = 'mov';
elseif ~isVis && ~isMov,
myTypes{ii} = 'none';
end
end

%% Plot these out
close all;
figure();
vInds = find(ismember(myTypes,'vis'));
for iv = 1:length(vInds),
    subplot(1,2,1);
    plot(allTimeCell{1}(allTimeCell{1} >= -300 & allTimeCell{1} <= 500),goodSDF{1}(vInds(iv),allTimeCell{1} >= -300 & allTimeCell{1} <= 500),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
    subplot(1,2,2);
    plot(allTimeCell{2}(allTimeCell{2} >= -500 & allTimeCell{2} <= 300),goodSDF{2}(vInds(iv),allTimeCell{2} >= -500 & allTimeCell{2} <= 300),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
end
% subplot(1,2,1);
% plot(allTimeCell{1}(allTimeCell{1} >= -300 & allTimeCell{1} <= 500),nanmean(goodSDF{1}(vInds(iv),allTimeCell{1} >= -300 & allTimeCell{1} <= 500),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
% subplot(1,2,2);
% plot(allTimeCell{2}(allTimeCell{2} >= -500 & allTimeCell{2} <= 300),goodSDF{2}(vInds(iv),allTimeCell{2} >= -500 & allTimeCell{2} <= 300),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
vSp1=subplot(1,2,1);
set(gca,'XLim',[-300,500],'XTick',-200:200:400,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','left');
vSp2=subplot(1,2,2);
set(gca,'XLim',[-500,300],'XTick',-400:200:200,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','right');
myYs = [get(vSp1,'YLim'); get(vSp2,'YLim')];
linkaxes([vSp1,vSp2],'y');
set(vSp1,'YLim',[min(myYs(:,1)),max(myYs(:,2))]);
vst=suptitle('Visual Unit SDFs'); set(vst,'fontsize',24);

figure();
mInds = find(ismember(myTypes,'mov'));
for iv = 1:length(mInds),
    subplot(1,2,1);
    plot(allTimeCell{1}(allTimeCell{1} >= -300 & allTimeCell{1} <= 500),goodSDF{1}(mInds(iv),allTimeCell{1} >= -300 & allTimeCell{1} <= 500),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
    subplot(1,2,2);
    plot(allTimeCell{2}(allTimeCell{2} >= -500 & allTimeCell{2} <= 300),goodSDF{2}(mInds(iv),allTimeCell{2} >= -500 & allTimeCell{2} <= 300),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
end
mSp1=subplot(1,2,1);
set(gca,'XLim',[-300,500],'XTick',-200:200:400,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','left');
mSp2=subplot(1,2,2);
set(gca,'XLim',[-500,300],'XTick',-400:200:200,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','right');
myYs = [get(mSp1,'YLim'); get(mSp2,'YLim')];
linkaxes([mSp1,mSp2],'y');
set(mSp1,'YLim',[min(myYs(:,1)),max(myYs(:,2))]);
mst=suptitle('Movement Unit SDFs'); set(mst,'fontsize',24);

figure();
vmInds = find(ismember(myTypes,'vismov'));
for iv = 1:length(vmInds),
    subplot(1,2,1);
    plot(allTimeCell{1}(allTimeCell{1} >= -300 & allTimeCell{1} <= 500),goodSDF{1}(vmInds(iv),allTimeCell{1} >= -300 & allTimeCell{1} <= 500),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
    subplot(1,2,2);
    plot(allTimeCell{2}(allTimeCell{2} >= -500 & allTimeCell{2} <= 300),goodSDF{2}(vmInds(iv),allTimeCell{2} >= -500 & allTimeCell{2} <= 300),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
end
vmSp1=subplot(1,2,1);
set(gca,'XLim',[-300,500],'XTick',-200:200:400,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','left');
vmSp2=subplot(1,2,2);
set(gca,'XLim',[-500,300],'XTick',-400:200:200,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','right');
myYs = [get(vmSp1,'YLim'); get(vmSp2,'YLim')];
linkaxes([vmSp1,vmSp2],'y');
set(vmSp1,'YLim',[min(myYs(:,1)),max(myYs(:,2))]);
vmst=suptitle('Visuomovement Unit SDFs'); set(vmst,'fontsize',24);

figure();
nInds = find(ismember(myTypes,'none'));
for iv = 1:length(nInds),
    subplot(1,2,1);
    plot(allTimeCell{1}(allTimeCell{1} >= -300 & allTimeCell{1} <= 500),goodSDF{1}(nInds(iv),allTimeCell{1} >= -300 & allTimeCell{1} <= 500),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
    subplot(1,2,2);
    plot(allTimeCell{2}(allTimeCell{2} >= -500 & allTimeCell{2} <= 300),goodSDF{2}(nInds(iv),allTimeCell{2} >= -500 & allTimeCell{2} <= 300),'color',[.2 .2 .2]+[.1 .1 .1].*mod(iv,3)); hold on;
end
nSp1=subplot(1,2,1);
set(gca,'XLim',[-300,500],'XTick',-200:200:400,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','left');
nSp2=subplot(1,2,2);
set(gca,'XLim',[-500,300],'XTick',-400:200:200,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','right');
myYs = [get(nSp1,'YLim'); get(nSp2,'YLim')];
linkaxes([nSp1,nSp2],'y');
set(nSp1,'YLim',[min(myYs(:,1)),max(myYs(:,2))]);
nst=suptitle('Unidentified Unit SDFs'); set(nst,'fontsize',24);

figure(); colors = 'rbgcm';
uTypes = unique(myTypes);
for iu = 1:length(uTypes),
    subplot(1,2,1);
%     pltMeanStd(allTimeCell{1}(allTimeCell{1} >= -300 & allTimeCell{1} <= 500),nanmean(goodSDF{1}(ismember(myTypes,uTypes{iu}),allTimeCell{1} >= -300 & allTimeCell{1} <= 500),1),nanstd(goodSDF{1}(ismember(myTypes,uTypes{iu}),allTimeCell{1} >= -300 & allTimeCell{1} <= 500),[],1)./sqrt(sum(ismember(myTypes,uTypes{iu}))),'color',colors(iu));
    plot(allTimeCell{1}(allTimeCell{1} >= -300 & allTimeCell{1} <= 500),nanmean(goodSDF{1}(ismember(myTypes,uTypes{iu}),allTimeCell{1} >= -300 & allTimeCell{1} <= 500),1),'color',colors(iu),'linewidth',2); hold on;
    subplot(1,2,2);
%     pltMeanStd(allTimeCell{2}(allTimeCell{2} >= -500 & allTimeCell{2} <= 300),nanmean(goodSDF{2}(ismember(myTypes,uTypes{iu}),allTimeCell{2} >= -500 & allTimeCell{2} <= 300),1),nanstd(goodSDF{2}(ismember(myTypes,uTypes{iu}),allTimeCell{2} >= -500 & allTimeCell{2} <= 300),[],1)./sqrt(sum(ismember(myTypes,uTypes{iu}))),'color',colors(iu));
    plot(allTimeCell{2}(allTimeCell{2} >= -500 & allTimeCell{2} <= 300),nanmean(goodSDF{2}(ismember(myTypes,uTypes{iu}),allTimeCell{2} >= -500 & allTimeCell{2} <= 300),1),'color',colors(iu),'linewidth',2); hold on;
end
aSp1=subplot(1,2,1);
set(gca,'XLim',[-300,500],'XTick',-200:200:400,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','left');
aSp2=subplot(1,2,2);
set(gca,'XLim',[-500,300],'XTick',-400:200:200,'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'box','off','YAxisLocation','right');
myYs = [get(aSp1,'YLim'); get(aSp2,'YLim')];
linkaxes([aSp1,aSp2],'y');
set(aSp1,'YLim',[min(myYs(:,1)),max(myYs(:,2))]);
ast=suptitle('Group Average SDFs'); set(ast,'fontsize',24);


%% Make confusion matrix with waveform IDs
% [p,x2,df,exp,obs] = klConting(myTypes,checkClustIDs);
% realDiff = (obs-exp)./exp;
% for i = 1:1000,
%     [p,x2,df,exp,obs] = klConting(klShuffle(myTypes),klShuffle(checkClustIDs));
%     shuffDiff(:,:,i) = (obs-exp)./exp;
% end
% mnDiff = nanmean(shuffDiff,3);
% stdDiff = nanstd(shuffDiff,[],3);
% zDiff = (realDiff-mnDiff)./stdDiff;
% figure();
% imagesc(zDiff);
% cb = colorbar; caxis([-3,3]);
% colormap(customMap);
% set(gca,'box','off','XTick',1:5,'YTick',1:4,'YTickLabel',unique(myTypes),'tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',18,'XTickLabel','');
% set(cb,'Ticks',-3:1.5:3,'tickdir','out','ticklength',get(cb,'ticklength').*3);


