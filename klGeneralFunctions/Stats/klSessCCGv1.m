function [allCCG, allTimes] = klSessCCGv1(spikesA,spikesB)

allCCG = [];
for ir = 1:size(spikesA,1)
    [ccg,ccgTimes,zInd] = klCCGv1(spikesA(ir,:),spikesB(ir,:));
    if isempty(allCCG),
        allCCG = ccg; currZ = zInd;
    else
        if currZ < zInd,
            allCCG = cat(2,nan(size(allCCG,1),zInd-currZ),allCCG,nan(size(allCCG,1),zInd-currZ));
            currZ = zInd;
        else
            ccg = cat(2,nan(1,currZ-zInd),ccg,nan(1,currZ-zInd));
        end
%         keyboard
        allCCG = cat(1,allCCG,ccg);
    end
end
allTimes = (-currZ+1):(currZ-1);