function saccs = klSaccDetect(inX,inY,varargin)

% Set defaults
minDur = 20;
smooth = 3;
nSD = 1;

% Decode varargin
varStrInd = 


if ~exist('thresh','var'),
    thresh = nSD*nanstd(diffX);
end

% Smooth, if asked:
if smooth,
    kern = klMakeGauss(smooth); kern = kern./sum(kern);
    smoothX = conv2(inX,kern,'same');
    smoothY = conv2(inY,kern,'same');
else
    smoothX = inX;
    smoothY = inY;
end

% Get derivatives to check for monotonic changes
diffX = diff(smoothX);
diffY = diff(smoothY);

% Check for increases (positive diff*)
% We want to find the time when its continuously increasing (posX) but only
% when a given threshold had been reached (posXThresh)
[posX, posXstarts] = klGetConsecutive(diffX > 0);
[posXThresh] = klGetConsecutive(diffX > nSD*nanstd(diffX));
posXSacc = unique(posXstarts(posXThresh > minDur));
posXEnd = posX(posXSacc)+posXSacc;

% Same for positive Y
[posY, posYstarts] = klGetConsecutive(diffY > nSD*nanstd(diffY));
[posYThresh] = klGetConsecutive(diffY > nSD*nanstd(diffY));
posYSacc = unique(posYstarts(posYThresh > minDur));
posYEnd = posY(posYSacc)+posYSacc;

% Now negative
[negX, negXstarts] = klGetConsecutive(diffX < 0);
[negXThresh] = klGetConsecutive(diffX < -(nSD*nanstd(diffX)));
negXSacc = unique(negXstarts(negXThresh > minDur));
negXEnd = negX(negXSacc)+negXSacc;

[negY, negYstarts] = klGetConsecutive(diffY < 0);
[negYThresh] = klGetConsecutive(diffY < -(nSD*nanstd(diffY)));
negYSacc = unique(negYstarts(negYThresh > minDur));
negYEnd = negY(negYSacc)+negYSacc;

% Now, we'll pull out the times where there were sustained
% increases/decreases for longer than minDur samples
allX = [posXSacc,negXSacc]; allXEnd = [posXEnd,negXEnd];
[xSacc,xInds] = sort(allX); xEnds = allXEnd(xInds);
allY = [posYSacc,negYSacc]; allYEnd = [posYEnd,negYEnd];
[ySacc,yInds] = sort(allY); yEnds = allYEnd(yInds);

% Now eliminate doubles
diffMat = abs(bsxfun(@minus,xSacc,ySacc'));
if any(diffMat(:) < minDur),
    [yDub,xDub] = find(diffMat < minDur);
    yDubVals = ySacc(yDub);
    xDubVals = xSacc(xDub);
    yDubVals(yDubVals <= xDubVals) = nan;
    xDubVals(xDubVals < yDubVals) = nan;
    yEnds(ismember(ySacc,yDubVals)) = [];
    xEnds(ismember(xSacc,xDubVals)) = [];
    ySacc(ismember(ySacc,yDubVals)) = [];
    xSacc(ismember(xSacc,xDubVals)) = [];
end

saccs = nan(length(ySacc) + length(xSacc),2);
allSaccs = [xSacc,ySacc]; allEnds = [xEnds,yEnds];
[saccs(:,1), sortInds] = sort(allSaccs);
saccs(:,2) = allEnds(sortInds);
