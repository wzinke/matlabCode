function klPlotClusts(tmpK,clustIDs,visSDF,visTimes,movSDF,movTimes)

close all;
for i = 1:tmpK,
    figure(i);
    subplot(1,2,1);
    plot(visTimes,visSDF(clustIDs(:,tmpK)==i,:),'color',[.7 .3 .3]);
    hold on;
    plot(visTimes,nanmean(visSDF(clustIDs(:,tmpK)==i,:),1),'color','k','linewidth',2);
    set(gca,'XLim',[-300,500]);
    
    subplot(1,2,2);
    plot(movTimes,movSDF(clustIDs(:,tmpK)==i,:),'color',[.7 .3 .3]);
    hold on;
    plot(movTimes,nanmean(movSDF(clustIDs(:,tmpK)==i,:),1),'color','k','linewidth',2);
    set(gca,'XLim',[-500,300]);
end