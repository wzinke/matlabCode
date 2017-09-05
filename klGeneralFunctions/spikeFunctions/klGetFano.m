function fano = klGetFano(spiketimes,varargin)

spkCounts = sum(isfinite(spiketimes),2);
fano = (nanstd(spkCounts).^2)./nanmean(spkCounts);