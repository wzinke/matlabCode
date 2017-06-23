function [outSorts] = tdtSortChanv2a(chanData,varargin)

warning off

% Set defaults
maxWvs = 10000;
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
[spkTimes,spkWaves,spkThresh] = klThreshCrossv4(chanData,'times',chanTimes,'-m',minRefract,'-t',thisThresh);

% Get a subset of the waveforms
if length(spkTimes) <= maxWvs,
    subInds = 1:length(spkTimes);
else
    randVect = randperm(maxWvs);
    subInds = randVect(1:maxWvs);
    
end
subWaves = spkWaves(subInds,:);
subTimes = spkTimes(subInds);
    
% Sort the waveforms
% [tmpSorts, alWaves, alTimes, scores, coeffs, eigVect] = klSortWaves(subWaves);

% Get aligned version of all 
if smooth && align,
    splWaves = spline(1:32,spkWaves,1:.1:32);
    [alWaves,alTimes] = klTroughAlignv4(splWaves,spline(1:32,wvTimes,1:.1:32),0);
else
    alWaves = spkWaves;
    alTimes = wvTimes;
end
subAligned = alWaves(subInds,:);

for ic = 1:size(alWaves,2),
    goodCols(ic) = sum(isfinite(alWaves(:,ic)))==size(alWaves,1);
end

% Get the PCA/LPP scores and do sorts
[coeffs,subScoresPCA] = pca(subAligned(:,goodCols),'Centered',0);
[subScoresLPP,eigVect] = klLPPv1(subAligned(:,goodCols));

[outK(1), outIDs{1}, outGap(1,:), outErr(1,:)] = klAutoSortv3a(subScoresPCA(:,1:nDims),'-p',1,'-t','pca','perc',1);
[outK(2), outIDs{2}, outGap(2,:), outErr(2,:)] = klLSCwGap(subScoresPCA(:,1:nDims),'-p',min([ceil(size(subScoresPCA,1)/2),1000]),'perc',1);
[outK(3), outIDs{3}, outGap(3,:), outErr(3,:)] = klAutoSortv3a(subScoresLPP(:,1:nDims),'-p',1,'-t','pca','perc',1);
[outK(4), outIDs{4}, outGap(4,:), outErr(4,:)] = klLSCwGap(subScoresLPP(:,1:nDims),'-p',min([ceil(size(subScoresLPP,1)/2),1000]),'perc',1);

gapComp = [outGap(:,2:end)-outErr(:,2:end),nan(size(outGap,1),1)];
for ig = 1:4,
    threshGap(ig) = outGap(ig,outK(ig));
    threshComp(ig) = gapComp(ig,outK(ig));
end
threshDiff = threshGap-threshComp;

tmpSorts.allKs = outK;
tmpSorts.idx = outIDs;
tmpSorts.gap = outGap;
tmpSorts.err = outErr;
tmpSorts.set = find(threshDiff == max(threshDiff));
tmpSorts.spkNum = 1;
tmpSorts.wvTime = alTimes;
tmpSorts.k = tmpSorts.allKs(tmpSorts.set);
tmpSorts.spkAvg = nanmean(alWaves(outIDs{tmpSorts.set}(:,tmpSorts.k)==tmpSorts.spkNum,:),1);
tmpSorts.isSig = zeros(1,tmpSorts.k);

% Compare groups from subset to the whole set of waveforms
switch compType,
    % Get means of the groups
    case 'corr', % Assign by maximum correlation
        mnSorts = nan(tmpSorts.k,size(subAligned,2));
        for i = 1:tmpSorts.k,
            mnSorts(i,:) = nanmean(subAligned(tmpSorts.idx{tmpSorts.set}(:,tmpSorts.k)==i,goodCols),1);
        end
        allCorrs = corr(spkWaves',mnSorts');
        [~,outSorts.idx] = max(allCorrs,[],2);
    case 'euc', % Get minimum euclidean distance from centroid
        % Get transformation
        if ismember(tmpSorts.set,[1,2]),
            
            % If original data = score*coeff', then score =
            % data*inv(coeff')?
            allScores = alWaves*inv(coeffs');
            subScores = subScoresPCA;
        else
            allScores = alWaves*eigVect;
            subScores = subScoresLPP;
        end
        
        % Get centroids
        mnSorts = nan(tmpSorts.k,size(subScores,2));
        for i = 1:tmpSorts.k,
            mnSorts(i,:) = nanmean(subScores(tmpSorts.idx{tmpSorts.set}(:,tmpSorts.k)==i,:),1);
        end
        sortDist = EuDist2(allScores,mnSorts);
        
        % Get minimum centroid distance
        [~,outSorts.idx] = min(sortDist,[],2);
end      

for i = 1:length(unique(outSorts.idx)),
    outSorts.mnWaves(i,:) = nanmean(alWaves(outSorts.idx==i,:),1);
end
outSorts.subInds = subInds;
outSorts.subSorts = tmpSorts;
outSorts.lpp = subScoresLPP;
outSorts.pca = subScoresPCA;
outSorts.pcaCoeffs = coeffs;
outSorts.lppEig = eigVect;
outSorts.alTimes = alTimes;
outSorts.alWaves = alWaves;
outSorts.spkTimesAll = spkTimes;
outSorts.subTimes = subTimes;