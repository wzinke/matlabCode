function plotPAClusts(normResp,inTimes,sortIDs,respAdj,raw,myK)

close all;

myFig=figure();
for ip = 1:4
    ap(ip) = subplot(2,2,ip); hold on;
end

% for ip = 1:(length(inTimes)/2)
%     ap(ip) = subplot(2,(length(inTimes)/4),ip); hold on;
% end

xStr = {'Stimulus','Saccade','Reward'};
grpColors = jet(myK);
for ik = 1:myK
%     figure(myFig);
%     for ip = 1:(length(inTimes)/2)
%         axes(ap(ip));
%         plot(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1)-nanmean(normResp{ip+(length(inTimes)/2)}(sortIDs(:,myK)==ik,:),1),'linewidth',2,'color',grpColors(ik,:));
%         plot(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1)-nanmean(normResp{ip+(length(inTimes)/2)}(sortIDs(:,myK)==ik,:),1),'linewidth',2,'color',grpColors(ik,:));
%         
%         if ik == myK
%             xlabel(xStr{ip});
%             set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3);
%             if ip == 1
%                 ylabel('Normalized Firing Rate (Z)');
%             else
%                 
%                 set(gca,'YAxisLoc','right');
%             end
%         end
%     end
    axes(ap(1));
    plot(inTimes{1},nanmean(normResp{1}(sortIDs(:,myK)==ik,:),1)-nanmean(normResp{3}(sortIDs(:,myK)==ik,:),1),'color',grpColors(ik,:));
    axes(ap(2));
    plot(inTimes{2},nanmean(normResp{2}(sortIDs(:,myK)==ik,:),1)-nanmean(normResp{4}(sortIDs(:,myK)==ik,:),1),'color',grpColors(ik,:));
    axes(ap(3));
    plot(inTimes{5},nanmean(normResp{5}(sortIDs(:,myK)==ik,:),1)-nanmean(normResp{7}(sortIDs(:,myK)==ik,:),1),'color',grpColors(ik,:));
    axes(ap(4));
    plot(inTimes{6},nanmean(normResp{6}(sortIDs(:,myK)==ik,:),1)-nanmean(normResp{8}(sortIDs(:,myK)==ik,:),1),'color',grpColors(ik,:));
    drawnow;
    pause(.3);
    
    figure();
    sp(1) = subplot(2,2,1); hold on;
    plot(inTimes{1},nanmean(normResp{1}(sortIDs(:,myK)==ik,:),1),'color','k','linewidth',3);
    plot(inTimes{3},nanmean(normResp{3}(sortIDs(:,myK)==ik,:),1),'color','k','linewidth',1);
    ylabel('Pro SDF');
    title('Stim Aligned');
    sp(2) = subplot(2,2,2); hold on;
    plot(inTimes{2},nanmean(normResp{2}(sortIDs(:,myK)==ik,:),1),'color','k','linewidth',3);
    plot(inTimes{4},nanmean(normResp{4}(sortIDs(:,myK)==ik,:),1),'color','k','linewidth',1);
    title('Sacc Aligned');
    sp(3) = subplot(2,2,3); hold on;
    plot(inTimes{5},nanmean(normResp{5}(sortIDs(:,myK)==ik,:),1),'color','r','linewidth',3);
    plot(inTimes{7},nanmean(normResp{7}(sortIDs(:,myK)==ik,:),1),'color','r','linewidth',1);
    ylabel('Anti SDF');
    sp(4) = subplot(2,2,4); hold on;
    plot(inTimes{6},nanmean(normResp{6}(sortIDs(:,myK)==ik,:),1),'color','r','linewidth',3);
    plot(inTimes{8},nanmean(normResp{8}(sortIDs(:,myK)==ik,:),1),'color','r','linewidth',1);
    linkaxes(sp,'y');
    suptitle(sprintf('Cluster %d - n=%d',ik,sum(sortIDs(:,myK)==ik)));
    drawnow;
    
    pause(.3);
    
%     for ip = 1:(length(inTimes)/2)
%         sp(ip) = subplot(1,(length(inTimes)/2),ip); hold on;
%     end
%     for ip = 1:(length(inTimes)/2)
%         axes(sp(ip));
%         pltMeanStd(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1),nanstd(normResp{ip}(sortIDs(:,myK)==ik,:),[],1)./sqrt(sum(sortIDs(:,myK)==ik)),'color',[.8 .2 .2]);
%         pltMeanStd(inTimes{ip},nanmean(normResp{ip+(length(inTimes)/2)}(sortIDs(:,myK)==ik,:),1),nanstd(normResp{ip+(length(inTimes)/2)}(sortIDs(:,myK)==ik,:),[],1)./sqrt(sum(sortIDs(:,myK)==ik)),'color',[.2 .2 .8]);
% %         plot(inTimes{ip},nanmean(normResp{ip}(sortIDs(:,myK)==ik,:),1),'color',[.8 .2 .2],'linewidth',2);
% %         plot(inTimes{ip},nanmean(normResp{ip+2}(sortIDs(:,myK)==ik,:),1),'color',[.2 .2 .8],'linewidth',2);
%         xlabel(xStr{ip});
%         set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3);
%         if ip == 1
%             ylabel('Normalized Firing Rate (Z)');
%         elseif ip == (length(inTimes)/2)
%             set(gca,'YAxisLoc','right');
%         else
%             set(gca,'ycolor',[1 1 1]);
%         end
%     end
%     suptitle(sprintf('Cluster %d - n=%d',ik,sum(sortIDs(:,myK)==ik)));
%     linkaxes(sp,'y');
end
% linkaxes(ap,'y');
        
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