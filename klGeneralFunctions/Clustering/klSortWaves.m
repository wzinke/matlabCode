function [chanSorts, alWaves, alTimes, scores, coeffs, eigVect] = klSortWaves(inWaves,varargin)

% Set defaults
spkFreq = 1000/24414;
wvTimes = ((1:32)-9).*spkFreq;
doLive = 0;
nDims =2;

% Decode varargin

cutWaves = klCutForSplinesv1(inWaves);
clear splWaves alWaves alTimes
for iw = 1:size(cutWaves,1),
    splWaves(iw,:) = spline(1:32,cutWaves(iw,:),1:.1:32);
end
[alWaves,alTimes] = klTroughAlignv4(splWaves,spline(1:32,wvTimes,1:.1:32),0);

for ic = 1:size(alWaves,2),
    goodCols(ic) = sum(isfinite(alWaves(:,ic)))==size(alWaves,1);
end

[outK(1),outIDs{1},outGap(1,:),outErr(1,:),scores{1},coeffs{1},eigVect{1}] = klPlotSortCombos(alWaves(:,goodCols),'pca','kmeans',1,0,'ndims',nDims);
suptitle('PCA/KMeans');

[outK(2),outIDs{2},outGap(2,:),outErr(2,:),scores{2},coeffs{2},eigVect{2}] = klPlotSortCombos(alWaves(:,goodCols),'pca','lsc',2,0,'ndims',nDims);
suptitle('PCA/LSC');

[outK(3),outIDs{3},outGap(3,:),outErr(3,:),scores{3},coeffs{3},eigVect{3}] = klPlotSortCombos(alWaves(:,goodCols),'lpp','kmeans',3,0,'ndims',nDims);
suptitle('LPP/KMeans');

[outK(4),outIDs{4},outGap(4,:),outErr(4,:),scores{4},coeffs{4},eigVect{4}] = klPlotSortCombos(alWaves(:,goodCols),'lpp','lsc',4,0,'ndims',nDims);
suptitle('LPP/LSC');

gapComp = [outGap(:,2:end)-outErr(:,2:end),nan(size(outGap,1),1)];
for ig = 1:4,
    threshGap(ig) = outGap(ig,outK(ig));
    threshComp(ig) = gapComp(ig,outK(ig));
end
threshDiff = threshGap-threshComp;

chanSorts.allKs = outK;
chanSorts.idx = outIDs;
chanSorts.gap = outGap;
chanSorts.err = outErr;

set = find(threshDiff == max(threshDiff));
spk = 1;
k = outK(set);
sig = 1;

if doLive,
    fprintf('Select set by typing "set = ..."\n');
    fprintf('Select spike group by typing "spk = ..."\n');
    fprintf('Change k by typing "k = ..."\n');

    keyboard
end

chanSorts.set = set;
chanSorts.spkNum = spk;
chanSorts.wvTime = alTimes;
chanSorts.k = k;
chanSorts.spkAvg = nanmean(alWaves(outIDs{set}(:,k)==spk,:),1);
chanSorts.isSig = sig;
