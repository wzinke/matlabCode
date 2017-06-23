function klPlotClustsv2(myClustIDs,visAlign, movAlign, seFlag)

global visTimes movTimes normWaves wvTimes
% myK = wvK;
% % myClustIDs = visClustIDs(:,respK == myK);
% myClustIDs = waveClustIDs(:,k == wvK);
doStats= 0;
plotIndiv = 0;
doSDF = 0;

uK = unique(myClustIDs); uK(isnan(uK)) = [];
myK = length(uK);

colors = jet(myK);
figure();
statSampRate = 5;
vTest = -200:statSampRate:300;
mTest = -300:statSampRate:200;

vSig = nan(1,range(vTest)/statSampRate);
mSig = nan(1,range(mTest)/statSampRate);
if doStats
    for i = 1:length(vSig),
        vSig(i) = anovan(visAlign(:,visTimes == vTest(i)),{myClustIDs},'display','off');
        mSig(i) = anovan(movAlign(:,movTimes == mTest(i)),{myClustIDs},'display','off');
    end
    vSig = pAdjust(vSig); mSig = pAdjust(mSig);
end

if doSDF
for i = 1:myK,
    subplot(1,2,1);
    if seFlag
        pltMeanStd(visTimes(ismember(visTimes,-200:300)),nanmean(visAlign(myClustIDs == uK(i),ismember(visTimes,-200:300)),1),nanstd(visAlign(myClustIDs == uK(i),ismember(visTimes,-200:300)),[],1)./sqrt(sum(myClustIDs == uK(i))),'color',colors(i,:));
    else
        plot(visTimes(ismember(visTimes,-200:300)),nanmean(visAlign(myClustIDs == uK(i),ismember(visTimes,-200:300)),1),'color',colors(i,:));
    end
    hold on; vline(0);
    if i == myK,
    %     ylabel('Normalized SDF','fontsize',18,'Rotation',0);
        xlabel('Time from StimOnset','fontsize',18);
        ylabel('Normalized SDF','fontsize',18);
        yVals=get(gca,'YLim');
        set(gca,'color','k');
        if doStats,
            scatter(vTest(vSig < .05),ones(1,sum(vSig < .05)).*(yVals(2)+(.1*range(yVals))),'r*');
        end
        lines(1) = hline(0); lines(2) = vline(0); set(lines,'color','w');
    end
    subplot(1,2,2);
    if seFlag,
        pltMeanStd(movTimes(ismember(movTimes,-300:200)),nanmean(movAlign(myClustIDs == uK(i),ismember(movTimes,-300:200)),1),nanstd(movAlign(myClustIDs == uK(i),ismember(movTimes,-200:300)),[],1)./sqrt(sum(myClustIDs == uK(i))),'color',colors(i,:));
    else
        plot(movTimes(ismember(movTimes,-300:200)),nanmean(movAlign(myClustIDs == uK(i),ismember(movTimes,-300:200)),1),'color',colors(i,:));
    end
    hold on; vline(0); hline(0);

    if i == myK,
        yVals=get(gca,'YLim');
        set(gca,'color','k');
        scatter(mTest(mSig < .05),ones(1,sum(mSig < .05)).*(yVals(2)+(.1*range(yVals))),'r*');
        xlabel('Time from SRT','fontsize',18);
        st=suptitle('Cluster Responses'); set(st,'fontsize',22);
        lines(1) = hline(0); lines(2) = vline(0); set(lines,'color','w');
    end
end
end

for i = 1:myK,
    if plotIndiv,
        figure(i+10);
        plot(wvTimes, normWaves(myClustIDs == uK(i),:)); hold on;
        plot(wvTimes,nanmean(normWaves(myClustIDs == uK(i),:),1),'k','linewidth',3);
    end
    
    figure(100);
    plot(wvTimes,nanmean(normWaves(myClustIDs == uK(i),:),1),'color',colors(i,:),'linewidth',3); hold on;
    legStr{i} = sprintf('n=%d',sum(myClustIDs == uK(i)));
end
legend(legStr,'color','w');
set(gca,'color','k');