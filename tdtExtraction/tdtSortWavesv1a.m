function [outSorts] = tdtSortWavesv1a(spkTimes,spkWaves,varargin)

warning off

% Set defaults
maxWvs = 40000;
maxAgglom = 2500;
minRefract = .6;
% thisThresh = -4*(median(abs(chanData)/.6745));
spkFreq = 1000/24414;
wvTimes = ((1:32)-9).*spkFreq;
compType = 'euc';
doPos = 0;
smooth = 0;
align = 0;
nDims = 2;
dimRedType = 'pca';
doTime = 1;
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
        case {'-d','nDims'},
            nDims = varargin{varStrInd(iv)+1};
        case {'-w'},
            maxWvs = varargin{varStrInd(iv)+1};
    end
end

% Get a subset of the waveforms
%if length(spkTimes) <= maxWvs,
%    subInds = 1:length(spkTimes);
%else
%    randVect = randperm(length(spkTimes));
%    subInds = randVect(1:maxWvs);
%end
if length(spkTimes) <= 30000,
    subInds = 1:length(spkTimes);
else
    randVect = randperm(length(spkTimes));
    subInds = randVect(1:30000);
end
subWaves = spkWaves(subInds,:);
subTimes = spkTimes(subInds);
    
% Get aligned version of all 
if align,
    [alWaves,alTimes] = klTroughAlignv4(subWaves,wvTimes,0);
else
    alWaves = subWaves;
    alTimes = wvTimes;
end

%% First deal with negative threshold crossings
for ic = 1:size(alWaves,2),
    goodCols(ic) = sum(isfinite(alWaves(:,ic)))==size(alWaves,1);
end

% Get the PCA/LPP scores and do sorts
switch dimRedType,
    case {'pca'},
        [coeffs,subScores] = pca(alWaves(:,goodCols),'Centered',0);
        outSorts.reConstruct = coeffs;
    case {'lpp'},
        [subScores,eigVect] = klLPPv1(alWaves(:,goodCols));
        outSorts.reConstruct = eigVect;
end

normObs = (subScores(:,1:nDims)-repmat(nanmean(subScores(:,1:nDims),1),size(subScores,1),1))./repmat(nanstd(subScores(:,1:nDims),[],1),size(subScores,1),1);
distros = ones(1,nDims);
if doTime,
    distros = cat(2,distros,2);
    normObs = cat(2,normObs,subTimes'./(max(subTimes)./nanmean(range(normObs,1))));
end


%% Do the clustering
[outK, idx, ~, ~, ~, randStruct] = klAgglomClust_newv4(normObs,'-distro',distros);
ksig = zeros(1,outK);
for ik = 1:outK,
    ksig(ik) = klGetSNRv1(alWaves(idx(:,outK)==ik,:)) >= minSNR;
end

outSorts.subInds = subInds;
outSorts.idx = idx;
%outSorts.aggMap = aggMap;
outSorts.allWaves = spkWaves;
outSorts.k = outK;
outSorts.isSig = ksig;
outSorts.subScores = subScores;
outSorts.dimRed = dimRedType;
outSorts.allTimes = spkTimes;
outSorts.randStruct = randStruct;