function [spiketimes, LFP, spkWaves] = tdtExtractSessv3(thisFile,events,rawDir,procDir),

if ~exist('rawDir','var') || isempty(rawDir),
    rawDir = 'Y:/Users/Kaleb/dataRaw';
end
if ~exist('procDir','var') || isempty(procDir),
    procDir = 'Y:/Users/Kaleb/dataProcessed';
end

% Set constants
lfpWind = [-500, 6500];

Task = tdtGetTaskv3(sprintf('%s/%s',rawDir,thisFile),events);
nTrs = length(Task.AlignTimes);
uLocs = unique(Task.TargetLoc(~isnan(Task.TargetLoc)));
colors = jet(length(uLocs));

% Save behavior only
fprintf('Saving...')
save(sprintf('%s/%s/Behav.mat',procDir,thisFile),'Task');
for ib = 1:length(['Saving...']), fprintf('\b'); end
fprintf('Done!\n');

% Now, channel by channel, get spikes and LFPs
numChans = length(dir(sprintf('%s/%s/*.sev',rawDir,thisFile)))/2;

fprintf('\t\tNeural Data: Found %d channels\n',numChans);
for ic = 1:numChans,
    fprintf('\t\t\tReading channel %d: ',ic);

%         clear chanSEVs maxSpk spiketimes LFP vSDF mSDF vTimes mTimes tmpFront tmpBack
    clearvars -except ic numChans colors uLocs nTrs Task procDir thisFile toDo events lfpWind lfpFreq spkFreq rawDir procDir

    % Read the data
    [spikes, LFP] = tdtGetChanv1(sprintf('%s/%s',rawDir,thisFile),ic);
    spiketimes = spikes.spiketimes;
    spkWaves = spikes.spkWaves;
    spkThresh = spikes.spkThresh;
    
    %% Let's make some summary plots
    % Visual responses
    figure(1); set(1,'visible','off');
    for ib = 1:length(['Getting LFP aligned on saccade...']), fprintf('\b'); end
    fprintf('Plotting visual response...');

    % Get SDF
    [vSDF,vTimes] = klSpkRatev2(spiketimes,'-q',1);

    % Plot spikes
    subplot(1,2,1);
    for il = 1:length(uLocs),
        plot(vTimes,nanmean(vSDF(Task.Correct==1 & Task.TargetLoc==uLocs(il),:),1),'color',colors(il,:)); hold on;
    end
    set(gca,'XLim',[-200 500]);
    xlabel('Time (ms)'); ylabel('Firing Rate (Hz)');

    % Plot LFP
    subplot(1,2,2);
    for il = 1:length(uLocs),
        plot(LFP.vTimes,nanmean(LFP.vis(Task.Correct==1 & Task.TargetLoc==uLocs(il),:),1),'color',colors(il,:)); hold on;
    end
    set(gca,'XLim',[-200 500]);
    xlabel('Time (ms)'); ylabel('LFP (?V)');

    suptitle(sprintf('%s - DSP%d - Visual',thisFile,ic));

    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s/DSP%d-visSummary.png',procDir,thisFile,ic));

    % Now movement responses
    figure(2); set(2,'visible','off');
    for ib = 1:length(['Plotting visual response...']), fprintf('\b'); end
    fprintf('Plotting saccade response...');

    % Get SDF
    [mSDF,mTimes] = klSpkRatev2(spiketimes-repmat(Task.GoCue+Task.SRT,1,size(spiketimes,2)),'-q',1);

    % Plot spikes
    subplot(1,2,1);
    for il = 1:length(uLocs),
        plot(mTimes,nanmean(mSDF(Task.Correct==1 & Task.TargetLoc==uLocs(il),:),1),'color',colors(il,:)); hold on;
    end
    set(gca,'XLim',[-500 200]);
    xlabel('Time (ms)'); ylabel('Firing Rate (Hz)');

    % Plot LFP
    subplot(1,2,2);
    for il = 1:length(uLocs),
        plot(LFP.mTimes,nanmean(LFP.mov(Task.Correct==1 & Task.TargetLoc==uLocs(il),:),1),'color',colors(il,:)); hold on;
    end
    set(gca,'XLim',[-500 200]);
    xlabel('Time (ms)'); ylabel('LFP (?V)');

    suptitle(sprintf('%s - Channel %d - Movement',thisFile,ic));

    % Save it
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s/DSP%d-movSummary.png',procDir,thisFile,ic));

    close(1); close(2);

    %% Save file data separately
    for ib = 1:length(['Plotting saccade response...']), fprintf('\b'); end
    fprintf('Saving channel...');

    % This channel spikes
    mkdir(sprintf('%s/%s/Channel%d',procDir,thisFile,ic));
    save(sprintf('%s/%s/Channel%d/Spikes.mat',procDir,thisFile,ic),'Task','spiketimes');

    % This channel LFPs
    save(sprintf('%s/%s/Channel%d/LFPs.mat',procDir,thisFile,ic),'Task','LFP');

    % This channel waveforms
    save(sprintf('%s/%s/Channel%d/Waves.mat',procDir,thisFile,ic),'spkWaves','spkThresh');

    for ib = 1:length(['Saving channel...']), fprintf('\b'); end
    fprintf('Done!\n');
    
    if ic < numChans,
        for ib = 1:length(sprintf('\t\t\tReading channel %d: ',ic)), fprintf('\b'); end
    end
end
