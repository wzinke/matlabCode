function [outSorts] = tdtSortChanv4a(chanData,varargin)

warning off

% Set defaults
maxWvs = 10000;
maxAgglom = 2500;
minRefract = .6;
% thisThresh = -4*(median(abs(chanData)/.6745));
thisThresh = -2*rms(chanData);
spkFreq = 1000/24414;
wvTimes = ((1:32)-9).*spkFreq;
chanTimes = 0:spkFreq:(spkFreq*(size(chanData,2)-1)); % 1000x multiplier converts to ms
compType = 'euc';
doPos = 0;
smooth = 0;
align = 0;
nDims = 2;
dimRedType = 'pca';
doTime = 0;
distros = ones(1,nDims);
minSNR = 1.25;

% Decode varargin
varStrInd=find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t'},
            thisThresh = varargin{varStrInd(iv)+1};
        case {'-p'},
            doPos = varargin{varStrInd(iv)+1};
        case {'-c'},
            compType = varargin{varStrInd(iv)+1};
    end
end

% Set threshold, get crossing times
[spkTimesN,spkWavesN,spkThreshN] = klThreshCrossv5(chanData,'times',chanTimes,'-m',minRefract,'-t',-thisThresh);
[spkTimesP,spkWavesP,spkThreshP] = klThreshCrossv5(chanData,'times',chanTimes,'-m',minRefract,'-t',thisThresh,'-p',1);
% clear chanData
chanData = [];

% Get a subset of the waveforms
if length(spkTimesN) <= maxWvs,
    subIndsN = 1:length(spkTimes);
else
    randVect = randperm(maxWvs);
    subIndsN = randVect(1:maxWvs);
end
subWavesN = spkWavesN(subIndsN,:);
subTimesN = spkTimesN(subIndsN);
    
% Get a subset of the waveforms
if length(spkTimesP) <= maxWvs,
    subIndsP = 1:length(spkTimesP);
else
    randVect = randperm(maxWvs);
    subIndsP = randVect(1:maxWvs);
end
subWavesP = spkWavesP(subIndsP,:);
subTimesP = spkTimesP(subIndsP);
 
% Sort the waveforms
% [outSorts, alWaves, alTimes, scores, coeffs, eigVect] = klSortWaves(subWaves);

% Get aligned version of all 
if align,
    [alWavesN,alTimesN] = klTroughAlignv4(subWavesN,wvTimes,0);
    [alWavesP,alTimesP] = klTroughAlignv4(subWavesP,wvTimes,0);
else
    alWavesN = subWavesN; alWavesP = subWavesP;
    alTimesN = wvTimes; alTimesP = wvTimes;
end

%% First deal with negative threshold crossings
for ic = 1:size(alWavesN,2),
    goodCols(ic) = sum(isfinite(alWavesN(:,ic)))==size(alWavesN,1);
end

% Get the PCA/LPP scores and do sorts
switch dimRedType,
    case {'pca'},
        [coeffs,subScores] = pca(alWavesN(:,goodCols),'Centered',0);
        outSorts.neg.reConstruct = coeffs;
    case {'lpp'},
        [subScores,eigVect] = klLPPv1(alWavesN(:,goodCols));
        outSorts.neg.reConstruct = eigVect;
end

normObs = (subScores(:,1:nDims)-repmat(nanmean(subScores(:,1:nDims),1),size(subScores,1),1))./repmat(nanstd(subScores(:,1:nDims),[],1),size(subScores,1),1);
distros = ones(1,nDims);
if doTime,
    distros = cat(2,distros,2);
    normObs = cat(2,normObs,subTimesN'./(max(subTimesN)./nanmean(range(normObs,1))));
end


%% Do the clustering
[outK, idx, aggMap, ~, ~, randStruct] = klAgglomClust_newv4(normObs,'-distro',distros);
ksig = zeros(1,outK);
for ik = 1:outK,
    ksig(ik) = klGetSNRv1(alWavesN(idx==ik,:)) >= minSNR;
end

outSorts.neg.subInds = subIndsN;
outSorts.neg.idx = idx;
outSorts.neg.aggMap = aggMap;
outSorts.neg.allWaves = spkWavesN;
outSorts.neg.k = outK;
outSorts.neg.isSig = ksig;
outSorts.neg.subScores = subScores;
outSorts.neg.dimRed = dimRedType;
outSorts.neg.allTimes = spkTimesN;
outSorts.neg.randStruct = randStruct;

%% Now do positive threshold crossings
clear goodCols coeffs subScores normObs outK idx aggMap randStruct
% [goodCols, coeffs, subScores, normObs, outK, idx, aggMap, randStruct] = deal([]);
for ic = 1:size(alWavesP,2),
    goodCols(ic) = sum(isfinite(alWavesP(:,ic)))==size(alWavesP,1);
end

% Get the PCA/LPP scores and do sorts
switch dimRedType,
    case {'pca'},
        [coeffs,subScores] = pca(alWavesP(:,goodCols),'Centered',0);
        outSorts.pos.reConstruct = coeffs;
    case {'lpp'},
        [subScores,eigVect] = klLPPv1(alWavesP(:,goodCols));
        outSorts.pos.reConstruct = eigVect;
end

normObs = (subScores(:,1:nDims)-repmat(nanmean(subScores(:,1:nDims),1),size(subScores,1),1))./repmat(nanstd(subScores(:,1:nDims),[],1),size(subScores,1),1);
if doTime,
    distros = cat(2,distros,2);
    normObs = cat(2,normObs,subTimesP'./(max(subTimesP)./nanmean(range(normObs,1))));
end

%% Do the clustering
[outK, idx, aggMap, ~, ~, randStruct] = klAgglomClust_newv4(normObs,'-distro',distros);
ksig = zeros(1,outK);
for ik = 1:outK,
    ksig(ik) = klGetSNRv1(alWavesP(idx==ik,:)) >= minSNR;
end

outSorts.pos.subInds = subIndsP;
outSorts.pos.idx = idx;
outSorts.pos.aggMap = aggMap;
outSorts.pos.allWaves = spkWavesP;
outSorts.pos.k = outK;
outSorts.pos.isSig = ksig;
outSorts.pos.subScores = subScores;
outSorts.pos.dimRed = dimRedType;
outSorts.pos.allTimes = spkTimesP;
outSorts.pos.randStruct = randStruct;
