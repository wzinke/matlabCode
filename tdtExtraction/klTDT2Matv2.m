function [Task, spiketimes, LFP, lfpOutTimes] = klTDT2Matv2(filepath,varargin),

%% Deal with options
% Set defaults
hasLFP = 1;
hasSpk = 1;
lfpWind = [-500, 6500];

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-e','events','ev'},
            events = varargin{varStrInd(iv)+1};
        case {'-c','chan'},
            doChans = varargin{varStrInd(iv)+1};
    end
end

% Assume events are from 028 unless otherwise specified
if ~exist('events','var'),
    events = TEMPO_EV_cosman_rig028;
end

%% Get the TDT structure using their TDT2mat
% First events
tdtEvsRaw = TDT2mat(filepath,'TYPE',{'epocs'});
tdtEvs = tdt2EvShft(tdtEvsRaw.epocs.EVNT.data);
tdtEvTms = tdtEvsRaw.epocs.EVNT.onset.*1000; % Multiplication converts to ms

% Now EEG
% tdtEEG = TDT2mat(filepath,'TYPE',{'streams'},'STORE','EEGG');

% SEV's (LFP and Wav1):
tdtSEVs = SEV2mat(filepath);

if isfield(tdtSEVs,'Lfp1'),
    tdtLFP = tdtSEVs.Lfp1.data;
    lfpRate = (24000/tdtSEVs.Lfp1.fs);
    lfpTimes = 0:lfpRate:(lfpRate*(size(tdtLFP,2)-1));
else
    fprintf('\tNo LFPs found on this session...\n');
    hasLFP = 0;
end

if isfield(tdtSEVs,'Wav1'),
    tdtSpk = tdtSEVs.Wav1.data;
    spkTimesRaw = 0:(1000/tdtSEVs.Wav1.fs):((1000/tdtSEVs.Wav1.fs)*(size(tdtSpk,2)-1)); % 1000x multiplier converts to ms
else
    fprintf('\tNo Spiking channels found on this session...\n');
    hasSpk = 0;
end

% Set channl  number
if ~exist('doChans','var') && hasSpk,
    doChans = 1:size(tdtSpk,1);
elseif ~exist('doChans','var')  && hasLFP,
    doChans = 1:size(tdtLFP,1);
end

% Let's set a threshold for spikes, if sorts are not included

% keyboard

% % %% Now let's start getting events and trials
% % evTypes = events.names;
% % for ie = 1:length(evTypes),
% %     evInds{ie} = find(tdtEvs == events.(evTypes{ie}));
% % end

% Let's get the number of trials:
trStartInds = find(tdtEvs == events.TrialStart_);
trStopInds  = find(tdtEvs == events.Eot_);
nTrs = length(trStartInds);

% We should get the maximum trial length for analog data...
trStarts = tdtEvTms(trStartInds);
trEnds   = [tdtEvTms(trStartInds(2:end));length(tdtEvTms)];
trTimes  = trEnds-trStarts;

% Get lengths of trials (in samples) for each trial
% maxLFP = ceil(max(trTimes)*tdtSEVs.Lfp1.fs);
% maxSpk = ceil(max(trTimes)*tdtSEVs.Wav1.fs);

% Let's loop through trials (can we make this a logical???)
% outLFPs = cell(size(tdtLFP,1));
% outSpks = cell(size(tdtSpk,1));
% for ic = 1:size(tdtLFP,1),
%     outLFPs{ic} = nan(nTrs,maxLFP);
%     outSpks{ic} = nan(nTrs,maxSpk);
% end

for it = 1:nTrs,
    % Let's see which code indices are good for this trial:
    strtInd = trStartInds(it);
    if it < nTrs,
        stopInd = trStartInds(it+1);
    else
        stopInd = length(tdtEvs);
    end
    thisTrInds = strtInd:stopInd;
%     thisTrLFPSamps = ceil(trTimes(it)*tdtSEVs.Lfp1.fs);
%     thisTrSpkSamps = ceil(trTimes(it)*tdtSEVs.Wav1.fs);
    
% % %     % Now let's loop through events:
% % %     for ie = 1:length(evTypes),
% % %         clear tmpInd
% % %         if any(ismember(evInds{ie},thisTrInds)),
% % %             tmpInd = thisTrInds(find(ismember(thisTrInds,evInds{ie}),1));
% % %             Task.(evTypes{ie})(it) = tdtEvTms(tmpInd);
% % %         else
% % %             Task.(evTypes{ie})(it) = nan;
% % %         end
% % %     end
% % %     
    % Now get LFP's and Spikes for this trial
%     for ic = 1:size(tdtLFP,1)
%         outLFPs{ic}(it,1:thisTrLFPSamps) = tdtLFP(ic,
%     end
end


%% Loop through trials to get trial info and make a trial info matrix (as per PLX conversion...)
% First get appropriate indices
cutTrs = [];
startInds = find(tdtEvs == events.StartInfos_)+1;
endInds = find(tdtEvs == events.EndInfos_)-1;
for it = 1:nTrs,
    
    if it < nTrs,
        tmpStart = find(startInds > trStopInds(it) & startInds < trStartInds(it+1),1,'first');
        tmpEnd   = find(endInds > trStopInds(it) & endInds < trStartInds(it+1),1,'first');
    else
        tmpStart = find(startInds > trStopInds(it),1,'first');
        tmpEnd = find(endInds > trStopInds(it),1,'first');
    end
    if isempty(tmpStart) || isempty(tmpEnd),
        cutTrs = [cutTrs;it];
    else
        iStart(it) = startInds(tmpStart);
        iEnd(it) = endInds(tmpEnd);
    end
    clear tmpStart tmpEnd
end
infoCnt = iEnd-iStart + 1;
evntCnt = trStopInds-trStartInds + 1;

% Now make the matrix
infoMat=nan(nTrs,max(infoCnt));
trialCodes = nan(nTrs,max(evntCnt));
trialTimes = nan(nTrs,max(evntCnt));
for it = 1:nTrs,
    infoMat(it,1:infoCnt(it)) = tdtEvs(iStart(it):iEnd(it));
    trialCodes(it,1:evntCnt(it)) = tdtEvs(trStartInds(it):trStopInds(it));
    trialTimes(it,1:evntCnt(it)) = tdtEvTms(trStartInds(it):trStopInds(it));
end

%% loop over trials to get relevant timings.

% Get photo diode event as stimulus onset, keep this time relative to trial
% onset for later time adjustments.
Task = decodeEvents(trialCodes,trialTimes,events,infoMat);
    
%% Now get spike events
if hasSpk,
    for ic = doChans,
        % Get sorts or set a threshold
        if exist('sortedChans','var'),
            spkTimes{ic} = sortedChans{ic};
        else
            [spkTimes{ic},spkWaves,spkThresh(ic)] = klThreshCrossv3(tdtSpk(ic,:),'times',spkTimesRaw);
        end

        % Count up spikes in each trial
        for it = 1:nTrs,
            numSpks(it,ic) = sum(spkTimes{ic} >= trStarts(it) & spkTimes{ic} <= trEnds(it));
        end
    end

    maxSpk = max(numSpks(:));
    spiketimes = nan(nTrs,maxSpk,max(doChans));

    % loop through channels and trials to place spikes
    for ic = doChans,
        for it = 1:nTrs,
            spiketimes(it,1:numSpks(it,ic),ic) =  spkTimes{ic}(spkTimes{ic} >= trStarts(it) & spkTimes{ic} <= trEnds(it));
        end
    end
    spiketimes = spiketimes - repmat(Task.AlignTimes,1,size(spiketimes,2),size(spiketimes,3));
else
    spiketimes = nan;
end

%% Finally, place LFPs
if hasLFP,
    LFP.vis = nan(nTrs,(length(lfpWind(1):lfpRate:lfpWind(2)))+1,max(doChans));
    for ic = doChans,
        for it = 1:nTrs,
            LFP.vis(it,1:sum(lfpTimes >= (Task.AlignTimes(it)+lfpWind(1)) & lfpTimes <= (Task.AlignTimes(it)+lfpWind(2))),ic) = tdtLFP(ic,lfpTimes >= (Task.AlignTimes(it)+lfpWind(1)) & lfpTimes <= (Task.AlignTimes(it)+lfpWind(2)));
        end
    end
    LFP.vTimes = ((1:size(LFP.vis,2))+lfpWind(1)).*lfpRate;
    
    % Now get it shifted and aligned on SRT
    for it = 1:nTrs,
        shiftInds(it) = max([0,find(LFP.vTimes >= (Task.SRT(it)+Task.GoCue(it)),1)]);
    end
    for ic = doChans
        tmpCell = mat2cell(LFP.vis(:,:,ic),ones(size(LFP.vis,1),1),size(LFP.vis,2));
        [LFP.mov(:,:,ic), newZero] = klAlignv5(tmpCell,shiftInds');
    end
    tmpFront = -1.*fliplr((0:1:newZero).*lfpRate);
    tmpBack = (1:(size(LFP.mov,2)-newZero-1)).*lfpRate;
    LFP.mTimes = [tmpFront,tmpBack];
    
else
    LFP = nan;
    lfpOutTimes = nan;
end
 
end
