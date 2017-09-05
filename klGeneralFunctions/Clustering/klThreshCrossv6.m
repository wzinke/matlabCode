function [outTimes, outMat, sortPols, thresh] = klThreshCrossv6(inVect,varargin)

% Set defaults
nSD = 1.25;
sampPoints = -8:23;
doPos = 0;
times = 1:length(inVect);
minRefract = 0;
checkPlus = 2;

% Decode varargin
varStrInd=find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','thresh'},
            thresh = varargin{varStrInd(iv)+1};
        case {'-s','sig'},
            nSD = varargin{varStrInd(iv)+1};
        case {'-p','pos'},
            doPos = varargin{varStrInd(iv)+1};
        case {'times'},
            times = varargin{varStrInd(iv)+1};
        case {'-m','-r'},
            minRefract = varargin{varStrInd(iv)+1};
    end
end

if ~exist('thresh','var'),
%     thresh = nSD*(nanstd(inVect));
    thresh = 4*(median(abs(inVect)/.6745)); % From Nguyen et al JNSciMeth 2014
end

sdLim = nSD*nanstd(inVect);

% Get positive threshold crossings
shift1 = [nan,inVect(1:(end-1))];
posCrosses = find(inVect > thresh & shift1 <= thresh);
posTimes = times(posCrosses);

% Get negative threshold crossings
negCrosses = find(inVect < -thresh & shift1 >= -thresh);
negTimes = times(negCrosses);

% Concatenate them
allCrosses = [posCrosses,negCrosses];
allTimes = [posTimes,negTimes];
allPols = [ones(1,length(posCrosses)),ones(1,length(negCrosses)).*2];

% Sort the times
[sortTimes,sortInds] = sort(allTimes);
sortCrosses = allCrosses(sortInds);
sortPols = allPols(sortInds);

% Cut final spike if it doesn't finish in time
sortPols(sortCrosses >= (length(inVect)-length(sampPoints))) = [];
sortCrosses(sortCrosses >= (length(inVect)-length(sampPoints))) = [];

% Now cut spikes that start too soon
sortPols(sortCrosses < (abs(sampPoints(1))+1)) = [];
sortCrosses(sortCrosses < (abs(sampPoints(1))+1)) = [];

if isempty(sortCrosses),
    outTimes = [];
    outMat = [];
    sortPols = [];
    return
end

outTimes = times(sortCrosses);
shiftTimes = [outTimes(2:end),nan];

% Now let's cut out ones where outTimes < minRefract
dTimes = [nan,diff(outTimes)];
%% get values "checkPlus" points later to see which to pick
% Get the indices where this is true
testInds = bsxfun(@plus,find(dTimes < minRefract)',[-1,0]);
% testInds = bsxfun(@plus,find(dTimes < minRefract)',[0,1]);
testVals(:,1) = max(abs(inVect(bsxfun(@plus,sortCrosses(testInds(:,1))',0:2))),[],2);
testVals(:,2) = max(abs(inVect(bsxfun(@plus,sortCrosses(testInds(:,2))',0:2))),[],2);

% Take the first if they are equal or both < nSD*sd(inVect)
sortCrosses(testInds(testVals(:,1)==testVals(:,2),2)) = nan;
sortCrosses(testInds(testVals(:,1) < sdLim & testVals(:,2) < sdLim,2)) = nan;

% Now, if one, but not the other, is below nSD*sd(inVect), keep the one
% that exceeds it
sortCrosses(testInds(testVals(:,1) < testVals(:,2) & testVals(:,1) < sdLim,1)) = nan;
sortCrosses(testInds(testVals(:,1) >= testVals(:,2) & testVals(:,2) < sdLim,2)) = nan;

% Finally, if they both exceed nSD*sd(inVect), take the first
sortCrosses(testInds(testVals(:,1) >= sdLim & testVals(:,2) >= sdLim,2)) = nan;

% Cut them out by taking away nans
outTimes(isnan(sortCrosses)) = [];
sortPols(isnan(sortCrosses)) = [];
sortCrosses(isnan(sortCrosses)) = [];

% As much as I'd like to take out the following loop, this is the most
% straightforward way I could think to minimize looping but get all
% necessary points...
outMat = single(nan(length(sortCrosses),length(sampPoints)));
for i = 1:length(sampPoints),
    outMat(:,i) = inVect(sortCrosses+sampPoints(i));
end

