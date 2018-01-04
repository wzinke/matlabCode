


rePull = 0;
whichType = 3;
if rePull
    [paSDF,paTimes,paTypes] = klPullProAnti;
else
    load('paSDFs-171026.mat');
end

% myCrit = find(paTypes(:,whichType));
myCrit = find(isnan(sortIDs(:,4)));
for row = 1:length(myCrit)
    fprintf('Plotting row %d of %d...\n',row,length(myCrit));
    ir = myCrit(row);
    figure();
    sp(1) = subplot(2,2,1);
    plot(paTimes{ir,1},paSDF{ir,1}(1,:),'color','k','linewidth',3);
    hold on;
    plot(paTimes{ir,1},paSDF{ir,1}(2,:),'color','k','linewidth',1);
    title('Stim Aligned');
    ylabel('Pro Firing Rate');
    subplot(2,2,2);
    sp(2) = subplot(2,2,2);
    plot(paTimes{ir,2},paSDF{ir,2}(1,:),'color','k','linewidth',3);
    hold on;
    plot(paTimes{ir,2},paSDF{ir,2}(2,:),'color','k','linewidth',1);
    title('Sacc Aligned');
    sp(3) = subplot(2,2,3);
    plot(paTimes{ir,1},paSDF{ir,1}(3,:),'color','r','linewidth',3);
    hold on;
    plot(paTimes{ir,1},paSDF{ir,1}(4,:),'color','r','linewidth',1);
    ylabel('Anti Firing Rate');
    sp(4) = subplot(2,2,4);
    plot(paTimes{ir,2},paSDF{ir,2}(3,:),'color','r','linewidth',3);
    hold on;
    plot(paTimes{ir,2},paSDF{ir,2}(4,:),'color','r','linewidth',1);
    linkaxes(sp,'y');
    
    pause(.5);
    saveas(gcf,sprintf('./anomaly%d.png',row));
    close all;
end