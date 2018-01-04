
% reDoMetaClustering
rePull = 0;
reClust = 0;
indivN = 10;
metaN = 10;
randReps = 1;

if rePull
    [goodSDF,goodSDFTimes] = klPullAllSDFs();
end
if reClust
    [sortIDs,~,rawM,respSumStruct] = klDoMetaClustering(goodSDF,goodSDFTimes,'-n',indivN,'-mn',metaN,'-r',randReps);
end

zDist = nan(size(goodSDF{1},1),size(goodSDF{1},1),length(respSumStruct));
for i = 1:length(respSumStruct)
    zDist(:,:,i) = (respSumStruct(i).distMat-nanmean(respSumStruct(i).distMat(:)))./nanstd(respSumStruct(i).distMat(:));
end

[idx,idxDist,~,~,linkMat,linkInds] = klDistMatAgglomv2(nanmedian(zDist,3),1:30,'-n',metaN);
% [idx,idxDist,~,~,linkMat,linkInds] = klDistMatAgglomv2(nanmean(zDist,3),1:30,'-n',metaN);
while sum(isnan(idx(:,end))) == size(idx,1)
    idx = idx(:,1:(end-1));
end
