function plotSearchClusts(normResp,inTimes,sortIDs,respAdj,raw,myK)

close all;

myFig=figure();
for ip = 1:(length(inTimes)/2)
    ap(ip) = subplot(1,(length(inTimes)/2),ip); hold on;
end
tColors = [.8 .2 .2; .8 .2 .2; .2 .8 .2];
dColors = [.2 .2 .8; .2 .2 .8; .8 .2 .8];
plotWindows = {[-100,300],[-200,200],[-200,200]};
xStr = {'Stimulus','Saccade','Reward'};
visCheck = 50:150;
movCheck = -25:25;
grpColors = jet(myK);
for ik = 1:myK
    figure(myFig);
    for ip = 1:(length(inTimes)/2)
        axes(ap(ip));
%         plot(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1)-nanmean(normResp{ip+(length(inTimes)/2)}(sortIDs(:,myK)==ik,:),1),'linewidth',2,'color',grpColors(ik,:));
        plot(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1),'linewidth',2,'color',grpColors(ik,:));
        if ik == myK
            xlabel(xStr{ip});
            set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'XLim',plotWindows{ip});
            if ip == 1
                ylabel('Normalized Firing Rate (Z)');
            else
                
                set(gca,'YAxisLoc','right');
            end
        end
    end
    clear sp
    figure();
    for ip = 1:(length(inTimes)/2)
        sp(ip) = subplot(1,(length(inTimes)/2),ip); hold on;
    end
    for ip = 1:(length(inTimes)/2)
        axes(sp(ip));
        pltMeanStd(inTimes{ip},nanmedian(normResp{ip}(sortIDs(:,myK)==ik,:),1),nanstd(normResp{ip}(sortIDs(:,myK)==ik,:),[],1)./sqrt(sum(sortIDs(:,myK)==ik)),'color',tColors(ip,:));
%         pltMeanStd(inTimes{ip},nanmedian(normResp{ip}(sortIDs(:,myK)==ik,:),1),nanstd(normResp{ip}(sortIDs(:,myK)==ik,:),[],1)./sqrt(sum(sortIDs(:,myK)==ik)),'color',grpColors(ik,:));
        pltMeanStd(inTimes{ip},nanmedian(normResp{ip+(length(inTimes)/2)}(sortIDs(:,myK)==ik,:),1),nanstd(normResp{ip+(length(inTimes)/2)}(sortIDs(:,myK)==ik,:),[],1)./sqrt(sum(sortIDs(:,myK)==ik)),'color',dColors(ip,:));
%         plot(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1),'color',[.8 .2 .2],'linewidth',2);
%         plot(inTimes{ip},nanmean(normResp{ip+2}(sortIDs(:,myK)==ik,:),1),'color',[.2 .2 .8],'linewidth',2);
        xlabel(xStr{ip});
        set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'XLim',plotWindows{ip});
        if ip == 1
            ylabel('Normalized Firing Rate (Z)');
        elseif ip == (length(inTimes)/2)
            set(gca,'YAxisLoc','right');
        else
            set(gca,'ycolor',[1 1 1],'YTick',[]);
        end
        
    end
    suptitle(sprintf('Cluster %d - n=%d',ik,sum(sortIDs(:,myK)==ik)));
    linkaxes(sp,'y');
    figure(100+ik);
    clear sp
    for ip = 1:(length(inTimes)/2)
        for it = 1:1
            sp(ip,it) = subplot(1,(length(inTimes)/2),ip+((length(inTimes)/2)*(it-1))); hold on;
        end
    end
    for ip = 1:(length(inTimes)/2)
        axes(sp(ip,1));
        plot(inTimes{ip},normResp{ip}(sortIDs(:,myK)==ik,:),'color',tColors(ip,:));
        plot(inTimes{ip},normResp{ip+(length(inTimes)/2)}(sortIDs(:,myK)==ik,:),'color',dColors(ip,:));
        xlabel(xStr{ip});
        set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'XLim',plotWindows{ip});
        if ip == 1
            ylabel('Normalized Firing Rate (Z)');
        elseif ip == (length(inTimes)/2)
            set(gca,'YAxisLoc','right');
        else
            set(gca,'ycolor',[1 1 1],'YTick',[]);
        end
    end
    suptitle(sprintf('Cluster %d - n=%d',ik,sum(sortIDs(:,myK)==ik)));
    linkaxes(sp(:,1),'y');
    
    figure(200+ik); hold on;
    myVisTarg = nanmean(normResp{1}(sortIDs(:,myK)==ik,ismember(inTimes{1},visCheck)),2);
    myVisDist = nanmean(normResp{1+(length(inTimes)/2)}(sortIDs(:,myK)==ik,ismember(inTimes{1+(length(inTimes)/2)},visCheck)),2);
    myMovTarg = nanmean(normResp{2}(sortIDs(:,myK)==ik,ismember(inTimes{2},movCheck)),2);
    myMovDist = nanmean(normResp{2+(length(inTimes)/2)}(sortIDs(:,myK)==ik,ismember(inTimes{2+(length(inTimes)/2)},movCheck)),2);
    scatter(myVisDist,myVisTarg,[],'k','filled');
    scatter(myMovDist,myMovTarg,[],'k');
    myX = get(gca,'XLim'); myY = get(gca,'YLim');
    plot([min([myX(1),myY(1)]),max([myX(2),myY(2)])],[min([myX(1),myY(1)]),max([myX(2),myY(2)])],'color',[.3 .3 .3],'linestyle','--');
    set(gca,'XLim',[min([myX(1),myY(1)]),max([myX(2),myY(2)])],'YLim',[min([myX(1),myY(1)]),max([myX(2),myY(2)])]);
    xlabel('Target Out'); ylabel('Target In');
end
linkaxes(ap,'y');
        
figure();
[h,t,perm] = dendrogram(respAdj,0,'Orientation','left');
klDendroClustChange(h,respAdj,sortIDs(:,myK));
set(gca,'XTick',[],'YTick',[],'box','off');

figure();
imagesc(-raw(fliplr(perm),fliplr(perm)));
set(gca,'XTick',[],'YTick',[]);
cb = colorbar;
set(cb,'tickdir','out');
colormap('jet');
caxis([-2 2]);