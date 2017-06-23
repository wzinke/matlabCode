%% 15-7-21 klGetMeanRate is a simple function that simply calculates the mean firing rate of a cell given spiketimes input
% Should work but doesn't seem to return the same as the mean of an sdf?
%       SDF could be dragged down by 0's/tails?

function [mr, sdr] = klGetMeanRate(spiketimes,varargin)

% Set constants
tConv = 1000; % tConv = time conversion. Assumed to be 1000 to convert ms to Hz
blWind = [-300,-100];
blOnly = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)},
        case {'base','-b'},
            blWind = varargin{varStrInd(iv)+1};
    end
end

if blOnly,
    spiketimes(spiketimes < blWind(1)) = nan;
    spiketimes(spiketimes > blWind(2)) = nan;
end

trSpks = sum(isfinite(spiketimes),2);
trTime = range(spiketimes,2);

% Cut nans (most importantly from time - no spike trials should be 0 in trSpks)
nanInd = isnan(trSpks) | isnan(trTime);
cumSpks = sum(trSpks(~nanInd)); 
% cumTime = sum(trTime(~nanInd));
if blOnly,
    cumTime = length(trTime).*diff(blWind);
else
    cumTime = nansum(trTime);
end
mr = cumSpks/(cumTime/tConv);

if nargout > 1,
   sdr = nanstd(trSpks(~nanInd)./(cumTime/tConv));
end