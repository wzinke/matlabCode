function [alWaves,  alTimes] = klSplAlignv1(waves,times)

diffWaves = [nan(size(waves,1),5),diff(waves(:,5:20),[],2),nan(size(waves,1),length(21:32))];
cutWaves=waves;
cutWaves(abs(diffWaves) < 2 & abs(waves) > nanstd(nanmean(waves,1))*2.5) = nan;
smoothWaves = nan(size(waves,1),length(1:.1:32));
for ii = 1:size(waves,1),
    smoothWaves(ii,:) = spline(1:32,cutWaves(ii,:),1:.1:32);
end
smoothTimes = spline(1:32,times,1:.1:32);
alignWind = 10;


[alWaves, alTimes] = klTroughAlignv4(smoothWaves,smoothTimes,0,'-w',alignWind);