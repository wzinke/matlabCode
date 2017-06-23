function [Task, trStarts, trEnds] = tdtGetTaskv2(inFile,events)

%% Get the TDT structure using their TDT2mat
% First events
tdtEvsRaw = TDT2mat(inFile,'TYPE',{'epocs'},'VERBOSE',false);
tdtEvs = tdt2EvShft(tdtEvsRaw.epocs.TEVT.data);
tdtEvTms = tdtEvsRaw.epocs.TEVT.onset.*1000; % Multiplication converts to ms

% Cut out weird zeros
tdtEvTms(tdtEvs <= 0) = [];
tdtEvs(tdtEvs <= 0) = [];

% Cut out codes that seem to just mess things up...
cutCodes = [3063, 4232, 8872, 4200, 128, 3997, 4032, 8064, 3976, 4224, 1024, 1094, 2176, 2944, 3072, 3968, 3584, 7936, 3456, 3328, 4096, 8832, 8192];
% cutCodes = [2176, 2944, 3072, 3584, 7936, 3456, 3328, 8832, 8192];
% cutCodes = [];

% Cut out duplicate codes that shouldn't be doubles...
keepDoubles = [8888,8200,8100];
tdtEvsShft = [tdtEvs(2:end);nan];
tdtEvs(tdtEvs==tdtEvsShft & ~ismember(tdtEvs,keepDoubles)) = nan;
tdtEvTms(isnan(tdtEvs)) = [];
tdtEvs(isnan(tdtEvs)) = [];

%% Loop through trials to get trial info and make a trial info matrix (as per PLX conversion...)
% First get appropriate indices
cutTrs = [];
startInds = find(tdtEvs == events.StartInfos_)+1;
endInds = find(tdtEvs == events.EndInfos_)-1;
trStartInds = find(tdtEvs == events.TrialStart_);
trStopInds  = find(tdtEvs == events.Eot_);
trStarts = tdtEvTms(trStartInds);
trEnds   = [tdtEvTms(trStartInds(2:end));length(tdtEvTms)];
nTrs = length(trStartInds);

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
    clear thisTrInfos
    thisTrInfos = tdtEvs(iStart(it):iEnd(it));
    % Cut out 2944 and 3072 which seems to be messing up the orders...
    thisTrInfos = thisTrInfos(~ismember(thisTrInfos,cutCodes)); 
    while sum(ismember(thisTrInfos(10:12),[8888,8200,8100])) == 3,
        if length(thisTrInfos) > 12,
            thisTrInfos = [thisTrInfos(1:11);thisTrInfos(13,:)];
        else
            thisTrInfos(12) = nan;
        end
    end
    trShift = [nan;thisTrInfos(1:(end-1))];
    thisTrInfos(thisTrInfos==trShift & ~ismember(thisTrInfos,[8888,8200,8100])) = [];
    if ismember(thisTrInfos(11),(0:45:315)+5000) && thisTrInfos(10) == 8888,
        thisTrInfos = [thisTrInfos(1:10);nan;thisTrInfos(11:end)];
    end
    
%     infoMat(it,1:infoCnt(it)) = tdtEvs(iStart(it):iEnd(it));
    infoMat(it,1:length(thisTrInfos)) = thisTrInfos;
    trialCodes(it,1:evntCnt(it)) = tdtEvs(trStartInds(it):trStopInds(it));
    trialTimes(it,1:evntCnt(it)) = tdtEvTms(trStartInds(it):trStopInds(it));
end

%% loop over trials to get relevant timings.


% Get photo diode event as stimulus onset, keep this time relative to trial
% onset for later time adjustments.
Task = decodeEvents(trialCodes,trialTimes,events,infoMat);
Task.trStarts = trStarts;
Task.trEnds = trEnds;
