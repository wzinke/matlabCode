function [spiketimes, LFP, spkWaves] = tdtExtractSessv4(thisFile,events,rawDir,procDir),

fresh = 0;

if ~exist('rawDir','var') || isempty(rawDir),
    rawDir = 'Y:/Users/Kaleb/dataRaw';
end
if ~exist('procDir','var') || isempty(procDir),
    procDir = 'Y:/Users/Kaleb/dataProcessed';
end

% Set constants
lfpWind = [-500, 2500];

Task = tdtGetTaskv4(sprintf('%s/%s',rawDir,thisFile),events);
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
    if fresh || ~(exist(sprintf('%s/%s/Channel%d/Spikes.mat',procDir,thisFile,ic),'file') && exist(sprintf('%s/%s/Channel%d/LFPs.mat',procDir,thisFile,ic),'file')),
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
        
        % Get SDF
        [vSDF,vTimes] = klSpkRatev2(spiketimes,'-q',1);

        % Plot spikes for MG Task
        if any(strcmp(Task.TaskType,'MG')),
            for ib = 1:length(['Getting LFP aligned on saccade...']), fprintf('\b'); end
            fprintf('Plotting MG visual response...');

            figure(1); set(1,'visible','off');
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
            for ib = 1:length(['Plotting MG visual response...']), fprintf('\b'); end
            fprintf('Plotting MG saccade response...');

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
            for ib = 1:length(['Plotting MG saccade response...']), fprintf('\b'); end

        end

        if any(ismember(Task.TaskType,{'Search','Cap'})),
            fprintf('Plotting Search visual response...');

            % Get locations
            capLocs = unique(Task.TargetLoc(ismember(Task.TaskType,{'Search','Cap'})));
            capLocs(isnan(capLocs)) = [];
            capColors = [.8 .2 .2; .2 .8 .2; .2 .2 .8];
            capDirs = 0:45:315;
            capMap = [6,3,2,1,4,7,8,9];

            % Plot out spikes
            figure(1); set(1,'visible','off');
            for il = 1:length(capLocs),
                subplot(3,3,capMap(capDirs==capLocs(il)));
                % Get indices for task type
                capInds(1,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc == capLocs(il);
                capInds(2,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= capLocs(il) & Task.Singleton==1 & Task.DistLoc == capLocs(il);
                capInds(3,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= capLocs(il) & (Task.Singleton==0 | Task.DistLoc == capLocs(il));

                % Plot them out
                for it = 1:3,
                    plot(vTimes,nanmean(vSDF(capInds(it,:),:),1),'color',capColors(it,:)); hold on;
                end
                set(gca,'XLim',[-200 500]);
                xlabel('Time (ms)');
            end
            suplabel('Firing Rate (Hz)','y');

            figure(2); set(2,'visible','off');
            % Plot LFP
            for il = 1:length(capLocs),
                subplot(3,3,capMap(capDirs==capLocs(il)));
                % Get indices for task type
                capInds(1,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc == capLocs(il);
                capInds(2,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= capLocs(il) & Task.Singleton==1 & Task.DistLoc == capLocs(il);
                capInds(3,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= capLocs(il) & (Task.Singleton==0 | Task.DistLoc == capLocs(il));

                % Plot them out
                for it = 1:3,
                    plot(LFP.vTimes,nanmean(LFP.vis(capInds(it,:),:),1),'color',capColors(it,:)); hold on;
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
            for ib = 1:length(['Plotting saccade response...']), fprintf('\b'); end

            figure(3); set(3,'visible','off');
            for ib = 1:length(['Plotting Search visual response...']), fprintf('\b'); end
            fprintf('Plotting Search saccade response...');

            % Get SDF
            [mSDF,mTimes] = klSpkRatev2(spiketimes-repmat(Task.GoCue+Task.SRT,1,size(spiketimes,2)),'-q',1);

            % Plot out spikes
            for il = 1:length(capLocs),
                subplot(3,3,capMap(capDirs==capLocs(il)));
                % Get indices for task type
                capInds(1,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc == capLocs(il);
                capInds(2,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= capLocs(il) & Task.Singleton==1 & Task.DistLoc == capLocs(il);
                capInds(3,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= capLocs(il) & (Task.Singleton==0 | Task.DistLoc == capLocs(il));

                % Plot them out
                for it = 1:3,
                    plot(mTimes,nanmean(mSDF(capInds(it,:),:),1),'color',capColors(it,:)); hold on;
                end
                set(gca,'XLim',[-500 400]);
                xlabel('Time (ms)');
            end
            suplabel('Firing Rate (Hz)','y');

            figure(4); set(4,'visible','off');
            % Plot LFP
            for il = 1:length(capLocs),
                subplot(3,3,capMap(capDirs==capLocs(il)));
                % Get indices for task type
                capInds(1,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc == capLocs(il);
                capInds(2,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= capLocs(il) & Task.Singleton==1 & Task.DistLoc == capLocs(il);
                capInds(3,:) = ismember(Task.TaskType,{'Search','Cap'}) & Task.Correct == 1 & Task.TargetLoc ~= capLocs(il) & (Task.Singleton==0 | Task.DistLoc == capLocs(il));

                % Plot them out
                for it = 1:3,
                    plot(LFP.mTimes,nanmean(LFP.mov(capInds(it,:),:),1),'color',capColors(it,:)); hold on;
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
            for ib = 1:length(['Plotting Search saccade response...']), fprintf('\b'); end

        end

        %% Save file data separately
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
    else
        continue
    end
end
