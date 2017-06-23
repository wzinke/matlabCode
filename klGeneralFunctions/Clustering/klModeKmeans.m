function [outID, outDists, outInd] = klModeKmeans(inIDs,dists)

% Assume inIDs is an #obs x #reps matrix with identical k values
uID = unique(inIDs(:,1));

for ik = 1:length(uID),
    for ir = 1:size(inIDs,2),
        idCount(ik,ir) = sum(inIDs(:,ir) == uID(ik));
    end
end

% Sort inIDs
sortIDs = sort(idCount,1);

modVals = nan(1,length(uID));
movIDs  = sortIDs;
rawInds = ones(1,size(inIDs,2));
for ik = 1:length(uID),
    modVals(ik) = mode(movIDs(ik,:));
    movIDs = movIDs(:,movIDs(ik,:) == modVals(ik));
    rawInds = rawInds & (sortIDs(ik,:) == modVals(ik));
end

outInd = find(rawInds,1);
outID = inIDs(:,outInd);
outDists = dists(:,outInd);