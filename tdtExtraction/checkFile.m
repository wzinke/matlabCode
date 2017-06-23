function checkFile(inFile,breakOut)

close all;

% Open/Decode the data file
[Task, spiketimes, LFP] = klTDT2Matv2(inFile);
uLocs = unique(Task.TargetLoc);
colors = jet(length(uLocs));
    
% Get SDFs
for ic = 1:size(spiketimes,3),
    [vSDF,vTimes] = klSpkRatev2(spiketimes(:,:,ic));
    [mSDF,mTimes] = klSpkRatev2(spiketimes(:,:,ic)-repmat(Task.SRT+Task.GoCue,1,size(spiketimes,2)));

    if breakOut
        figure(ic); 
        % Open and plot graph of SDFs
        visAxSpk = axes('position',[.4 .1 .5 .8]);
        for i = 1:length(uLocs),
            plot(vTimes,nanmean(vSDF(Task.Correct==1 & Task.TargetLoc==uLocs(i),:),1),'color',colors(i,:)); hold on;
        end
        set(gca,'XLim',[-200 500]);

        % Put in corresponding LFPs
        visAxLFP = axes('position',[.1 .1 .2 .8]);
        for i = 1:length(uLocs),
            plot(LFP.vTimes,nanmean(LFP.vis(Task.Correct==1 & Task.TargetLoc==uLocs(i),:),1),'color',colors(i,:)); hold on;
        end
        set(gca,'XLim',[-200 500]);

        figure(ic+100);
        movAxSpk = axes('position',[.4 .1 .5 .8]);
        for i = 1:length(uLocs),
            plot(mTimes,nanmean(mSDF(Task.Correct==1 & Task.TargetLoc==uLocs(i),:),1),'color',colors(i,:)); hold on;
        end
        set(gca,'XLim',[-500 200]);

        % Get shifted LFPs
        movAxSpk = axes('position',[.1 .1 .2 .8]);
        for i = 1:length(uLocs),
            plot(LFP.mTimes,nanmean(LFP.mov(Task.Correct==1 & Task.TargetLoc==uLocs(i),:),1),'color',colors(i,:)); hold on;
        end
        set(gca,'XLim',[-500 200])

        if pause
            keyboard
        end
    end
    figure(1001);
%     visAxes{ic} = axes('position',[.1 .1 .2 .1]);
    subplot(8,4,ic);
    plot(vTimes,nanmean(vSDF(Task.Correct==1,:),1));
    set(gca,'XLim',[-200 500]);
    
    figure(1002);
    subplot(8,4,ic);
    plot(LFP.vTimes,nanmean(LFP.vis(Task.Correct==1,:,ic),1));
    set(gca,'XLim',[-200 500]);
    
    figure(1003);
    subplot(8,4,ic);
    plot(mTimes,nanmean(mSDF(Task.Correct==1,:),1));
    set(gca,'XLim',[-500 200]);
    
    figure(1004);
    subplot(8,4,ic);
    plot(LFP.mTimes,nanmean(LFP.mov(Task.Correct==1,:,ic),1));
    set(gca,'XLim',[-500 200]);
    
    
end
