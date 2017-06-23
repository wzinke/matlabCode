function allIDX = klExtrapSortv1(allScores,subScores,idx,simThresh)

if nargin < 4,
    simThresh = [];
end

checkMax = 5000;
changeMax = 0;
allIDX = uint8(zeros(size(allScores,1),1));
myK = length(unique(idx(idx < 31)));

for ic = 1:ceil(size(allScores,1)/checkMax),
    theseVals = (1:checkMax)+(checkMax*(ic-1)); theseVals(theseVals > size(allScores,1)) = [];
    checkDist = EuDist2(allScores(theseVals,:),subScores);
    [~,minInd] = min(checkDist,[],2);
    allIDX(theseVals) = uint8(idx(minInd));
end

% Now, assign "unsorted" ones if they are closer to a cluster than the
% similarity at which the combination is made
nChanged = inf;
while nChanged > changeMax,
    nChanged = sum(allIDX == 31);
    unSorted = find(allIDX==31);
    
    % Loop like above
    for ic = 1:ceil(nChanged/checkMax),
        checkKs = nan(myK,checkMax);
        theseVals = (1:checkMax)+(checkMax*(ic-1)); theseVals(theseVals > nChanged) = [];
        for ik = 1:myK,
            checkDist = EuDist2(allScores(allIDX==ik,:),allScores(unSort
