function [fVect, rVect, betVect, witVect] = klGetF(obs,groups,varargin)

% Set defaults
type = 'euc';
isSim = 0;
catVars = zeros(1,size(obs,2));
simMat = [];
doN = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','type'}
            type = varargin{varStrInd(iv)+1};
        case {'-s','sim'}
            simMat = varargin{varStrInd(iv)+1};
            isSim = 1;
        case {'cat'}
            catVars = varargin{varStrInd(iv)+1};
        case {'-n','n'},
            doN = varargin{varStrInd(iv)+1};
    end
end

obsRaw = obs;
if ismember(type,{'prod'}),
    maxObs = nanmax(obs,[],1);
    minObs = nanmin(obs,[],1);
    normObs = (obs-repmat(minObs,size(obs,1),1))./repmat((maxObs-minObs),size(obs,1),1);
    normRaw = normObs;
    obs = normObs;
end

% Make similarity matrix, if necessary
if isempty(simMat),
    fprintf('Computing pairwise distances...');
    for ir = 1:size(obs,1),
        simMat(ir,ir) = 0;
        for ic = ir:size(obs,1),
            thisDist = klGetDist(obs(ir,:),obs(ic,:),'-t',type,'cat',catVars);
            simMat(ir,ic) = thisDist;
            simMat(ic,ir) = thisDist;
        end
    end
    fprintf('Done!\n');
end

% Loop through each grouping segment
fVect = nan(1,size(groups,2));
rVect = nan(1,size(groups,2));
betVect = nan(1,size(groups,2));
witVect = nan(1,size(groups,2));
for ic = 1:size(groups,2),
    % Loop through each unique group to get within group similarities
    obsMean = nanmean(obs,1);
    uGroup = unique(groups(~isnan(groups(:,ic)),ic));
    sumWithin = zeros(1,length(uGroup));
    sumTot = zeros(1,length(uGroup));
    sumBetween = zeros(1,length(uGroup));
    for ig = 1:length(uGroup),
        thisGroup = uGroup(ig);
        grpMean = nanmean(obs(groups(:,ic) == thisGroup,:),1);
        myObs = find(groups(:,ic) == thisGroup);
        for in = 1:sum(groups(:,ic) == thisGroup),
            sumWithin(ig) = sumWithin(ig)+klGetDist(obs(myObs(in),~catVars),grpMean(~catVars),'-t',type);
            sumTot(ig)    = sumTot(ig)+klGetDist(obs(myObs(in),~catVars),obsMean(~catVars),'-t',type);
        end
        sumBetween(ig) = klGetDist(grpMean(~catVars),obsMean(~catVars),'-t',type);
        if doN,
            sumBetween(ig) = sumBetween(ig)*sum(groups(:,ic) == thisGroup);
        end
    end
    dfBetween = length(uGroup)-1;
    dfWithin  = size(obs,1)-length(uGroup);
    fVect(ic) = (sum(sumBetween)/dfBetween)/(sum(sumWithin)/dfWithin);
    rVect(ic) = sum(sumBetween)/sum(sumTot);
    betVect(ic) = sum(sumBetween);
    witVect(ic) = sum(sumWithin);
end