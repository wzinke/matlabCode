function [outSorts] = tdtSortChanv1a(chanData,varargin)

warning off

% Set defaults
maxWvs = 10000;
minRefract = .2;
thisThresh = -4*(median(abs(chanData)/.6745));
spkFreq = 1000/24414;
wvTimes = ((1:32)-9).*spkFreq;
chanTimes = 0:spkFreq:(spkFreq*(size(chanData,2)-1)); % 1000x multiplier converts to ms
compType = 'euc';
doPos = 0;
smooth = 0;
align = 0;

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
[spkTimes,spkWaves,spkThresh] = klThreshCrossv4(chanData,'times',chanTimes,'-p',doPos,'-m',minRefract,'-t',thisThresh);

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
tmpSorts = klSortWaves(subWaves);

% Get aligned version of all 
% cutWaves = klCutForSplinesv1(spkWaves);
% for iw = 1:size(cutWaves,1),
%     splWaves(iw,:) = spline(1:32,cutWaves(iw,:),1:.1:32);
% end
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

[coeffs,subScoresPCA] = pca(subAligned(:,goodCols),'Centered',0);
[subScoresLPP,eigVect] = klLPPv1(subAligned(:,goodCols));
                  
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
outSorts.alTimes = alTimes;
outSorts.alWaves = alWaves;
outSorts.spkTimesAll = spkTimes;
outSorts.subTimes = subTimes;