function plotMGClusts(normResp,inTimes,sortIDs,respAdj,raw,myK)

close all;

myFig=figure();
for ip = 1:length(normResp)
    ap(ip) = subplot(1,length(normResp),ip); hold on;
end
grpColors = jet(myK);

for ik = 1:myK
    figure(myFig);
    for ip = 1:length(normResp)
        axes(ap(ip));
        plot(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1),'color',grpColors(ik,:),'linewidth',2);
    end
    figure();
    for ip = 1:length(normResp)
        sp(ip) = subplot(1,length(normResp),ip); hold on;
    end
    for ip = 1:length(normResp)
        axes(sp(ip));
        plot(inTimes{ip},normResp{ip}(sortIDs(:,myK)==ik,:),'color',[.4 .4 .4]);
        plot(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1),'color','k','linewidth',3);
    end
    linkaxes(sp,'y');
    set(sp(1),'XLim',[-200,300]);
    set(sp(2),'XLim',[-300,200]);

end
linkaxes(ap,'y');
set(ap(1),'XLim',[-200,300]);
set(ap(2),'XLim',[-300,200]);

        
figure();
[h,t,perm] = dendrogram(respAdj,0,'Orientation','left');
klDendroClustChange(h,respAdj,sortIDs(:,myK));
set(gca,'XTick',[],'YTick',[],'box','off');

figure();
imagesc(-raw(fliplr(perm),fliplr(perm)));
cb = colorbar;
colormap('jet');
caxis([-2 2]);