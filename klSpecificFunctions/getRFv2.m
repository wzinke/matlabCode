function [rf, antirf] = getRFv2(spks,Task,varargin)

bl = 0;
vWind = [50,150];
mWind = [-100,0];
doSDF = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'bl','-b'}
            bl = varargin{varStrInd(iv)+1};
        case {'-s','sdf'},
            doSDF = varargin{varStrInd(iv)+1};
    end
end

tLocs = unique(Task.TargetLoc(isfinite(Task.TargetLoc)));
vfr = nan(1,length(tLocs));
mfr = nan(1,length(tLocs));
mspks = spks-repmat(Task.SRT+Task.GoCue,1,size(spks,2));
if doSDF,
    [spks,vtimes] = klSpkRatev2(spks,'-q',1);
    [mspks,mtimes] = klSpkRatev2(mspks,'-q',1);
    
    for il = 1:length(tLocs),
        vfr(il) = abs(nanmean(nanmean(spks(Task.TargetLoc == tLocs(il),vtimes >= vWind(1) & vtimes <= vWind(2)),1),2) - bl);
        mfr(il) = abs(nanmean(nanmean(mspks(Task.TargetLoc == tLocs(il),mtimes >= mWind(1) & mtimes <= mWind(2)),1),2) - bl);
    end
else
    for il = 1:length(tLocs)
        vfr(il) = abs(nanmean(sum(spks(Task.TargetLoc == tLocs(il),:) >= vWind(1) & spks(Task.TargetLoc == tLocs(il),:) <= vWind(2),2),1)-bl);
        mfr(il) = abs(nanmean(sum(mspks(Task.TargetLoc == tLocs(il),:) >= mWind(1) & mspks(Task.TargetLoc == tLocs(il),:) <= mWind(2),2),1)-bl);
    end
end
rfT{1}    = tLocs(vfr == max(vfr));
rfT{2}    = tLocs(mfr ==  max(mfr));
antirfT{1} = tLocs(vfr == min(vfr));
antirfT{2} = tLocs(mfr == min(mfr));

if ~isempty(rfT{1}),
    rf(1) = rfT{1}(1);
else
    rf(1) = nan;
end
if ~isempty(rfT{2}),
    rf(2) = rfT{2}(1);
else
    rf(2) = nan;
end
if ~isempty(antirfT{1}),
    antirf(1) = rfT{1}(1);
else
    antirf(1) = nan;
end
if ~isempty(antirfT{2}),
    antirf(2) = rfT{2}(1);
else
    antirf(2) = nan;
end
