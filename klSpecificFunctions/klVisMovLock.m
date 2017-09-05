function evLock = klVisMovLock(poissLats, srt)

minTr = 5;
alph = .05;

% Start by cutting out nans
nanTrs = isnan(poissLats) | isnan(srt);

% Break if there aren't enough trials
if sum(~nanTrs) <= minTr,
    evLock = nan;
    return;
end
goodLats = poissLats(~nanTrs);
goodSRT  = srt(~nanTrs);

% Correlate with SRT
[~,p] = corr(goodLats,goodSRT);
if p < alph, evLock = 'mov'; end

% Get event locking by variance
visVar = var(goodLats);
movVar = var(goodSRT-goodLats);

keyboard