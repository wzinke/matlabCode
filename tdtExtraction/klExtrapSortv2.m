function allIDX = klExtrapSortv2(allScores,subScores,idx,simThresh)

if nargin < 4,
    simThresh = [];
end

checkMax = 5000;
changeMax = 0;
allIDX = uint8(zeros(size(allScores,1),1));
myK = length(unique(idx(idx < 31)));
maxLoops = 500;

for ic = 1:ceil(size(allScores,1)/checkMax),
    theseVals = (1:checkMax)+(checkMax*(ic-1)); theseVals(theseVals > size(allScores,1)) = [];
    checkDist = EuDist2(allScores(theseVals,:),subScores);
    [~,minInd] = min(checkDist,[],2);
    allIDX(theseVals) = uint8(idx(minInd));
end

% Now, assign "unsorted" ones if they are closer to a cluster than the
% similarity at which the combination is made
nChanged = inf; nLoops = 0;
while nChanged > changeMax && nLoops <= maxLoops,
    unSorted = sum(allIDX == 31);
    unSortInd = find(allIDX==31);
    
    % Loop like above
    for ic = 1:ceil(unSorted/checkMax),
        theseVals = (1:checkMax)+(checkMax*(ic-1)); theseVals(theseVals > unSorted) = [];
        checkKs = nan(myK,length(theseVals));
        for ik = 1:myK,
            checkDist = EuDist2(allScores(allIDX==ik,:),allScores(unSortInd(theseVals),:));
            checkKs(ik,:) = min(checkDist,[],1);
        end
        % If this observation is closer to more than one group, don't
        % change it
        checkKs(:,(sum(checkKs < simThresh,1) > 1)) = deal(nan);
        [minVals,minInd] = min(checkKs,[],1);
        allIDX(unSortInd(theseVals(minVals < simThresh))) = uint8(minInd(minVals < simThresh));
        nChanged = sum(minVals < simThresh);
    end
    nLoops = nLoops+1;
end
