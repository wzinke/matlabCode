function [Task, trStarts, trEnds] = tdtGetTaskv1(inFile,events)

%% Get the TDT structure using their TDT2mat
% First events
tdtEvsRaw = TDT2mat(inFile,'TYPE',{'epocs'},'VERBOSE',false);
tdtEvs = tdt2EvShft(tdtEvsRaw.epocs.EVNT.data);
tdtEvTms = tdtEvsRaw.epocs.EVNT.onset.*1000; % Multiplication converts to ms

% Cut out weird zeros
tdtEvTms(tdtEvs==min(tdtEvs)) = [];
tdtEvs(tdtEvs==min(tdtEvs)) = [];


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
    infoMat(it,1:infoCnt(it)) = tdtEvs(iStart(it):iEnd(it));
    trialCodes(it,1:evntCnt(it)) = tdtEvs(trStartInds(it):trStopInds(it));
    trialTimes(it,1:evntCnt(it)) = tdtEvTms(trStartInds(it):trStopInds(it));
end

%% loop over trials to get relevant timings.

% Get photo diode event as stimulus onset, keep this time relative to trial
% onset for later time adjustments.
Task = decodeEvents(trialCodes,trialTimes,events,infoMat);
