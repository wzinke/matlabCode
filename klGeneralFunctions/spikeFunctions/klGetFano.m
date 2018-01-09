function fano = klGetFano(spiketimes,varargin)

bin = 100;
spkN = hist(spiketimes,0:bin:max(spiketimes(:)));
fano = nanvar(spkN(:))./nanmean(spkN(:));
% spkCounts = sum(isfinite(spiketimes),2);
% fano = (nanstd(spkCounts).^2)./nanmean(spkCounts);