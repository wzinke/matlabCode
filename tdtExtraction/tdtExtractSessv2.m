function [spiketimes, LFP, spkWaves] = tdtExtractSessv2(thisFile,events,rawDir,procDir),

if ~exist('rawDir','var') || isempty(rawDir),
    rawDir = 'Y:/Users/Kaleb/dataRaw';
end
if ~exist('procDir','var') || isempty(procDir),
    procDir = 'Y:/Users/Kaleb/dataProcessed';
end

% Set constants
lfpWind = [-500, 6500];
lfpFreq = 1000/1017; % 1000 converts to ms
spkFreq = 1000/24414;

[Task, trStarts, trEnds] = tdtGetTaskv2(sprintf('%s/%s',rawDir,thisFile),events);
nTrs = length(Task.AlignTimes);
uLocs = unique(Task.TargetLoc(~isnan(Task.TargetLoc)));
colors = jet(length(uLocs));

% Save behavior only
fprintf('Saving...')
save(sprintf('%s/%s/Behav.mat',procDir,thisFile),'Task','trStarts','trEnds');
for ib = 1:length(['Saving...']), fprintf('\b'); end
fprintf('Done!\n');

% Now, channel by channel, get spikes and LFPs
numChans = length(dir(sprintf('%s/%s/*.sev',rawDir,thisFile)))/2;

fprintf('\t\tNeural Data: Found %d channels\n',numChans);
for ic = 1:numChans,
    fprintf('\t\t\tReading channel %d: ',ic);

%         clear chanSEVs maxSpk spiketimes LFP vSDF mSDF vTimes mTimes tmpFront tmpBack
    clearvars -except ic numChans colors uLocs nTrs Task trStarts trEnds procDir thisFile toDo events lfpWind lfpFreq spkFreq rawDir procDir

    % Read the data
    chanSEVs = SEV2mat_kl(sprintf('%s/%s',rawDir,thisFile),'CHANNEL',ic,'VERBOSE',0);
    
    lfpRate = (24000/chanSEVs.Lfp1.fs);
    lfpTimes = 0:lfpRate:(lfpRate*(size(chanSEVs.Lfp1.data,2)-1));

    spkFreq = 1000/24414;

    spkTimesRaw = 0:spkFreq:(spkFreq*(size(chanSEVs.Wav1.data,2)-1)); % 1000x multiplier converts to ms
%     lfpTimes = 0:lfpFreq:(lfpFreq*(size(chanSEVs.Lfp1.data,2)-1));

     % Set a threshold or load sort data
    fprintf('Getting spike times...');
    if exist('sortedChans','var'),
        spkTimes = sortedChans{ic};
    else
        [spkTimes,spkWaves,spkThresh] = klThreshCrossv3(chanSEVs.Wav1.data,'times',spkTimesRaw);
    end

    for ib = 1:length(['Getting spike times...']), fprintf('\b'); end
    fprintf('Counting spikes...');
    % Figure out how big the matrix will be
    nSpks = nan(nTrs,1);
    for it = 1:nTrs,
        nSpks(it) = sum(spkTimes >= trStarts(it) & spkTimes <= trEnds(it));
    end
    maxSpk = max(nSpks(:));
    spiketimes = nan(nTrs,maxSpk);

    for ib = 1:length(['Counting spikes...']), fprintf('\b'); end
    fprintf('Placing spikes...');
    % Place spikes in the matrix
    for it = 1:nTrs,
        spiketimes(it,1:nSpks(it)) =  spkTimes(spkTimes >= trStarts(it) & spkTimes <= trEnds(it));
    end
    spiketimes = spiketimes-repmat(Task.AlignTimes,1,size(spiketimes,2));

    for ib = 1:length(['Placing spikes...']), fprintf('\b'); end
    fprintf('Getting visually aligned LFP...');

    % Now get LFP
    LFP.vis = nan(nTrs,(length(lfpWind(1):lfpRate:lfpWind(2)))+1);
%     for ic = doChans,
        for it = 1:nTrs,
            LFP.vis(it,1:sum(lfpTimes >= (Task.AlignTimes(it)+lfpWind(1)) & lfpTimes <= (Task.AlignTimes(it)+lfpWind(2)))) = chanSEVs.Lfp1.data(1,lfpTimes >= (Task.AlignTimes(it)+lfpWind(1)) & lfpTimes <= (Task.AlignTimes(it)+lfpWind(2)));
        end
%     end
    LFP.vTimes = ((1:size(LFP.vis,2))+lfpWind(1)).*lfpRate;
    
    % Now get it shifted and aligned on SRT
    for ib = 1:length(['Getting visually aligned LFP...']), fprintf('\b'); end
    fprintf('Getting LFP aligned on saccade...');

    for it = 1:nTrs,
        shiftInds(it) = max([0,find(LFP.vTimes >= (Task.SRT(it)+Task.GoCue(it)),1)]);
    end
    tmpCell = mat2cell(LFP.vis,ones(size(LFP.vis,1),1),size(LFP.vis,2));
    [LFP.mov, newZero] = klAlignv5(tmpCell,shiftInds');

    tmpFront = -1.*fliplr((0:1:newZero).*lfpFreq);
    tmpBack = (1:(size(LFP.mov,2)-newZero-1)).*lfpFreq;
    LFP.mTimes = [tmpFront,tmpBack];

    %% Let's make some summary plots
    % Visual responses
    figure(1); set(1,'visible','off');
    for ib = 1:length(['Getting LFP aligned on saccade...']), fprintf('\b'); end
    fprintf('Plotting visual response...');

    % Get SDF
    [vSDF,vTimes] = klSpkRatev2(spiketimes,'-q',1);

    % Plot spikes for MG Task
    if any(strcmp(Task.TaskType,'MG')),
        
        subplot(1,2,1);
        for il = 1:length(uLocs),
            plot(vTimes,nanmean(vSDF(strcmp(Task.TaskType,'MG') & Task.Correct==1 & Task.TargetLoc==uLocs(il),:),1),'color',colors(il,:)); hold on;
        end
        set(gca,'XLim',[-200 500]);
        xlabel('Time (ms)'); ylabel('Firing Rate (Hz)');

        % Plot LFP
        subplot(1,2,2);
        for il = 1:length(uLocs),
            plot(LFP.vTimes,nanmean(LFP.vis(strcmp(Task.TaskType,'MG') & Task.Correct==1 & Task.TargetLoc==uLocs(il),:),1),'color',colors(il,:)); hold on;
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
            plot(mTimes,nanmean(mSDF(strcmp(Task.TaskType,'MG') & Task.Correct==1 & Task.TargetLoc==uLocs(il),:),1),'color',colors(il,:)); hold on;
        end
        set(gca,'XLim',[-500 200]);
        xlabel('Time (ms)'); ylabel('Firing Rate (Hz)');

        % Plot LFP
        subplot(1,2,2);
        for il = 1:length(uLocs),
            plot(LFP.mTimes,nanmean(LFP.mov(strcmp(Task.TaskType,'MG') & Task.Correct==1 & Task.TargetLoc==uLocs(il),:),1),'color',colors(il,:)); hold on;
        end
        set(gca,'XLim',[-500 200]);
        xlabel('Time (ms)'); ylabel('LFP (?V)');

        suptitle(sprintf('%s - Channel %d - Movement',thisFile,ic));

        % Save it
        set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(gcf,sprintf('%s/%s/DSP%d-movSummary.png',procDir,thisFile,ic));

        close(1); close(2);
    end
    
    if any(ismember(Task.TaskType,{'Search,'Cap'})),
        % Get locations
        capLocs = unique(Task.TargetLoc(ismember(Task.TaskType,{'Search','Cap'})));
        capColors = 'rgb';
        capDirs = 0:45:315;
        
        % Plot out spikes
        figure(1);
        for il = 1:length(capLocs),
            subplot(3,3,find(capDirs==capLocs(il)));
            % Get indices for task type
            capInds(1,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc == uLocs(il);
            capInds(2,:) = ismember(Task.taskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= uLocs(il) & Task.Singleton==1;
            capInds(3,:) = ismember(Task.taskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= uLocs(il) & Task.Singleton==0;
            
            % Plot them out
            for it = 1:3,
                plot(vTimes,nanmean(vSDF(capInds,1),'color',colors(it,:)); hold on;
            end
            set(gca,'XLim',[-200 500]);
            xlabel('Time (ms)');
        end
        suplabel('Firing Rate (Hz)','y');
        
        figure(2);
        % Plot LFP
        for il = 1:length(capLocs),
            subplot(3,3,find(capDirs==capLocs(il)));
            % Get indices for task type
            capInds(1,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc == uLocs(il);
            capInds(2,:) = ismember(Task.taskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= uLocs(il) & Task.Singleton==1;
            capInds(3,:) = ismember(Task.taskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= uLocs(il) & Task.Singleton==0;
            
            % Plot them out
            for it = 1:3,
                plot(LFP.vTimes,nanmean(LFP.vis(capInds,1),'color',colors(it,:)); hold on;
            end
            set(gca,'XLim',[-200 500]);
            xlabel('Time (ms)');
        end
        suplabel('LFP (?V)','y');

        suptitle(sprintf('%s - Search - DSP%d - Visual',thisFile,ic));

        set(1,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(1,sprintf('%s/%s/Search-DSP%d-visSummary-spks.png',procDir,thisFile,ic));
        set(2,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(2,sprintf('%s/%s/Search-DSP%d-visSummary-LFP.png',procDir,thisFile,ic));

        % Now movement responses
        figure(2); set(2,'visible','off');
        for ib = 1:length(['Plotting visual response...']), fprintf('\b'); end
        fprintf('Plotting saccade response...');

        % Get SDF
        [mSDF,mTimes] = klSpkRatev2(spiketimes-repmat(Task.GoCue+Task.SRT,1,size(spiketimes,2)),'-q',1);

        % Plot out spikes
        figure(1);
        for il = 1:length(capLocs),
            subplot(3,3,find(capDirs==capLocs(il)));
            % Get indices for task type
            capInds(1,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc == uLocs(il);
            capInds(2,:) = ismember(Task.taskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= uLocs(il) & Task.Singleton==1;
            capInds(3,:) = ismember(Task.taskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= uLocs(il) & Task.Singleton==0;
            
            % Plot them out
            for it = 1:3,
                plot(mTimes,nanmean(mSDF(capInds,1),'color',colors(it,:)); hold on;
            end
            set(gca,'XLim',[-500 400]);
            xlabel('Time (ms)');
        end
        suplabel('Firing Rate (Hz)','y');
        
        figure(4);
        % Plot LFP
        for il = 1:length(capLocs),
            subplot(3,3,find(capDirs==capLocs(il)));
            % Get indices for task type
            capInds(1,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc == uLocs(il);
            capInds(2,:) = ismember(Task.taskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= uLocs(il) & Task.Singleton==1;
            capInds(3,:) = ismember(Task.taskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= uLocs(il) & Task.Singleton==0;
            
            % Plot them out
            for it = 1:3,
                plot(LFP.mTimes,nanmean(LFP.mov(capInds,1),'color',colors(it,:)); hold on;
            end
            set(gca,'XLim',[-500 200]);
            xlabel('Time (ms)');
        end
        suplabel('LFP (?V)','y');

        suptitle(sprintf('%s - Search - DSP%d - Movement',thisFile,ic));

        set(3,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(3,sprintf('%s/%s/Search-DSP%d-movSummary-spks.png',procDir,thisFile,ic));
        set(4,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(4,sprintf('%s/%s/Search-DSP%d-movSummary-LFP.png',procDir,thisFile,ic));
        
        close all
    end
    
    %% Save file data separately
    for ib = 1:length(['Plotting saccade response...']), fprintf('\b'); end
    fprintf('Saving channel...');

    % This channel spikes
    mkdir(sprintf('%s/%s/Channel%d',procDir,thisFile,ic));
    save(sprintf('%s/%s/Channel%d/Spikes.mat',procDir,thisFile,ic),'Task','trStarts','trEnds','spiketimes');

    % This channel LFPs
    save(sprintf('%s/%s/Channel%d/LFPs.mat',procDir,thisFile,ic),'Task','trStarts','trEnds','LFP');

    % This channel waveforms
    save(sprintf('%s/%s/Channel%d/Waves.mat',procDir,thisFile,ic),'spkWaves','spkThresh');

    for ib = 1:length(['Saving channel...']), fprintf('\b'); end
    fprintf('Done!\n');
    
    if ic < numChans,
        for ib = 1:length(sprintf('\t\t\tReading channel %d: ',ic)), fprintf('\b'); end
    end
end
