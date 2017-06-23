function [outSorts] = tdtSortChanv3a(chanData,varargin)

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
align = 1;
nDims = 2;
dimRedType = 'pca';
doTime = 1;

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
if doTime,
    normObs = cat(2,normObs,subTimesN'./(max(subTimesN)./nanmean(range(normObs,1))));
end


%% Loop through subInds in manageable chunks
indsOpen = ones(1,size(subScores,1));

for is = 1:ceil(size(subScores,1)./maxAgglom),
    startInd = 1+(maxAgglom*(is-1));
    endInd = min([maxAgglom+(maxAgglom*(is-1)),length(subIndsN)]);
    currInds = subIndsN(startInd:endInd);
    [outKN{is},aggMembsN{is},clustColsN{is},outGapN{is}] = klAgglomClust_newv3(normObs(currInds,:));
    outIDsN{is} = ones(length(currInds),sum(isfinite(clustColsN{is}))).*31;
    for ic = 1:sum(isfinite(clustColsN{is})),
        lens = cellfun(@length,aggMembsN{is}(:,clustColsN{is}(ic)));
        [~,sortInds] = sort(lens,'descend');
        for ik = 1:ic,
            outIDsN{is}(aggMembsN{is}{sortInds(ik),clustColsN{is}(ic)},ic) = deal(ik);
        end
    end

    keyboard
end

outSorts.neg.subInds = subIndsN;
outSorts.neg.idx = outIDsN;
outSorts.neg.gap = outGapN;
outSorts.neg.spkNum = 1;
outSorts.neg.wvTime = alTimesN;
outSorts.neg.allWaves = spkWavesN;
outSorts.neg.k = outKN;
outSorts.neg.isSig = zeros(1,outKN);
outSorts.neg.subScores = subScores;
outSorts.neg.dimRed = dimRedType;
outSorts.neg.allTimes = spkTimesN;

%% Now do positive threshold crossings
clear goodCols coeffs subScores normObs
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
    normObs = cat(2,normObs,subTimesP'./(max(subTimesP)./nanmean(range(normObs,1))));
end
[outKP,aggMembsP,clustColsP,outGapP] = klAgglomClust_newv2(normObs);
outIDsP = ones(size(subScores,1),sum(isfinite(clustColsP))).*31;
for ic = 1:sum(isfinite(clustColsP)),
    lens = cellfun(@length,aggMembsP(:,clustColsP(ic)));
    [~,sortInds] = sort(lens,'descend');
    for ik = 1:ic,
        outIDsN(aggMembsP{sortInds(ik),clustColsP(ic)},ic) = deal(ik);
    end
end

outSorts.pos.subInds = subIndsP;
outSorts.pos.idx = outIDsP;
outSorts.pos.gap = outGapP;
outSorts.pos.spkNum = 1;
outSorts.pos.wvTime = alTimesP;
outSorts.pos.allWaves = spkWavesP;
outSorts.pos.k = outKP;
outSorts.pos.isSig = zeros(1,outKP);
outSorts.pos.subScores = subScores;
outSorts.pos.dimRed = dimRedType;
outSorts.pos.allTimes = spkTimesP;