function rf = getRF(spks,Task)

tLocs = unique(Task.TargetLoc(isfinite(Task.TargetLoc)));
fr = nan(1,length(tLocs));
for il = 1:length(tLocs)
    fr(il) = nanmean(sum(isfinite(spks(Task.TargetLoc == tLocs(il),:)),2),1);
end

rf = tLocs(fr == max(fr));