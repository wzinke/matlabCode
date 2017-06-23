function klTDT2Mat(filepath,varargin),

%% Deal with options
% Set defaults
hasLFP = 1;
hasSpk = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-e','events','ev'},
            events = varargin{varStrInd(iv)+1};
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
tdtEvTms = tdtEvsRaw.epocs.EVNT.onset;

% Now EEG
% tdtEEG = TDT2mat(filepath,'TYPE',{'streams'},'STORE','EEGG');

% SEV's (LFP and Wav1):
tdtSEVs = SEV2mat(filepath);

if isfield(tdtSEVs,'Lfp1'),
    tdtLFP = tdtSEVs.Lfp1.data;
else
    fprintf('\tNo LFPs found on this session...\n');
    hasLFP = 0;
end

if isfield(tdtSEVs,'Wav1'),
    tdtSpk = tdtSEVs.Wav1.data;
else
    fprintf('\tNo Spiking channels found on this session...\n');
    hasSpk = 0;
end

% Let's set a threshold for spikes, if sorts are not included

% keyboard

%% Now let's start getting events and trials
evTypes = events.names;
for ie = 1:length(evTypes),
    evInds{ie} = find(tdtEvs == events.(evTypes{ie}));
end

% Let's get the number of trials:
trStartInds = find(tdtEvs == events.TrialStart_);
nTrs = length(trStartInds);

% We should get the maximum trial length for analog data...
trStarts = tdtEvTms(trStartInds);
trEnds   = [tdtEvTms(trStartInds(2:end));length(tdtEvTms)];
trTimes  = trEnds-trStarts;

% Get lengths of trials (in samples) for each trial
maxLFP = ceil(max(trTimes)*tdtSEVs.Lfp1.fs);
maxSpk = ceil(max(trTimes)*tdtSEVs.Wav1.fs);

% Let's loop through trials (can we make this a logical???)
outLFPs = cell(size(tdtLFP,1));
outSpks = cell(size(tdtSpk,1));
for ic = 1:size(tdtLFP,1),
    outLFPs{ic} = nan(nTrs,maxLFP);
    outSpks{ic} = nan(nTrs,maxSpk);
end

for it = 1:nTrs,
    % Let's see which code indices are good for this trial:
    strtInd = trStartInds(it);
    if it < nTrs,
        stopInd = trStartInds(it+1);
    else
        stopInd = length(tdtEvs);
    end
    thisTrInds = strtInd:stopInd;
    thisTrLFPSamps = ceil(trTimes(it)*tdtSEVs.Lfp1.fs);
    thisTrSpkSamps = ceil(trTimes(it)*tdtSEVs.Wav1.fs);
    
    % Now let's loop through events:
    for ie = 1:length(evTypes),
        clear tmpInd
        if any(ismember(evInds{ie}),thisTrInds),
            tmpInd = thisTrInds(find(ismember(thisTrInds,evInds{ie}),1));
            outEvs.(evTypes{ie})(it) = tdtEvTms(tmpInd);
        else
            outEvs.(evTypes{ie})(it) = nan;
        end
    end
    
    % Now get LFP's and Spikes for this trial
    for ic = 1:size(tdtLFP,1)
%         outLFPs{ic}(it,1:thisTrLFPSamps) = tdtLFP(ic,
    end
end


