function [starts,ends] = klIdx2Perm(idx,perm),

uIDs = unique(idx(~isnan(idx)));
for iu = 1:length(uIDs),
    myVals = find(ismember(perm,find(idx==uIDs(iu))));
    starts(iu) = min(myVals);
    ends(iu) = max(myVals);
end
