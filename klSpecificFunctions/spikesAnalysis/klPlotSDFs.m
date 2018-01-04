function klPlotSDFs(spikeMat,waves,Task)

% spikeMat = klPlaceEvents(Task,spikes);
mgTrials = ismember(Task.TaskType,'MG');


[vSDF,vTimes] = klSpkRatev2(spikeMat,'-q',1);
[mSDF,mTimes] = klSpkRatev2(spikeMat-repmat(Task.SRT+Task.GoCue,1,size(spikeMat,2)),'-q',1);

close all;

figure();
sp(1) = subplot(3,2,[1,3]);
pltMeanStd(vTimes,nanmean(vSDF(mgTrials,:),1),nanstd(vSDF(mgTrials,:),[],1)./sqrt(sum(mgTrials)),'color','k');
set(gca,'XLim',[-300,500]);
sp(2) = subplot(3,2,[2,4]);
pltMeanStd(mTimes,nanmean(mSDF(mgTrials,:),1),nanstd(mSDF(mgTrials,:),[],1)./sqrt(sum(mgTrials)),'color','k');
set(gca,'XLim',[-500,300]);
linkaxes(sp,'y');

wp(1) = subplot(3,2,5);
r=randperm(size(waves,1));
plot(waves(r(1:min([length(r),10000])),:)');

wp(2) = subplot(3,2,6);
pltMeanStd(1:size(waves,2),nanmean(waves,1),nanstd(waves,[],1)./sqrt(size(waves,1)),'color','k');
linkaxes(wp,'y');
