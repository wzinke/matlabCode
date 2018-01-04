function [placedStreams, outTimes] = klPlaceStream(Task,stream,varargin)

% Set defaults
load('slowSamplingRate.mat');
pullRange = [-500, 3000];

timeVect = (0:(length(stream)-1)).*(1000/sampRate);

trStarts = Task.trStarts;
nTrs = length(trStarts);

startInd = find(timeVect > abs(pullRange(1)),1); startInd = startInd*((1i^(pullRange(1) < 0))^2);
endInd = find(timeVect < abs(pullRange(2)),1,'last'); endInd = endInd*((1i^(pullRange(2) < 0))^2);

placedStreams = nan(nTrs,length(startInd:endInd));
for it = 1:nTrs
    tmp = find(timeVect >= Task.AlignTimes(it),1);
    if ~isempty(tmp)
        placedStreams(it,:) = stream(tmp+(startInd:endInd)); 
    end
end
outTimes = (startInd:endInd).*(1000/sampRate);