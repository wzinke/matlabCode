function placedEvents = klPlaceEvents(Task,events)

maxSpkTm = 3000;


trStarts = Task.trStarts;
trEnds = Task.trEnds;
nTrs = length(trStarts);

nSpks = nan(nTrs,1);
for it = 1:nTrs
    nSpks(it) = sum(events >= trStarts(it) & events <= trEnds(it));
end
maxSpk = max(nSpks(:));
placedEvents = nan(nTrs,maxSpk);

% Place spikes in the matrix
for it = 1:nTrs,
    placedEvents(it,1:nSpks(it)) =  events(events >= trStarts(it) & events <= trEnds(it));
end
placedEvents = placedEvents-repmat(Task.AlignTimes,1,size(placedEvents,2));
placedEvents(placedEvents >= maxSpkTm) = nan;
cutCols = zeros(1,size(placedEvents,2));
for iCol = 1:size(placedEvents,2)
    cutCols(iCol) = sum(isnan(placedEvents(:,iCol)))==size(placedEvents,1);
end
placedEvents(:,logical(cutCols)) = [];