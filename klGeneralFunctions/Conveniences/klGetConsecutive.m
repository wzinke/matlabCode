function [out,starts] = klGetConsecutive(in)

out = zeros(size(in));
starts = nan(size(in));
for ir = 1:size(in,1),
    inRow = in(ir,:);
    nConsecRow = out(ir,:);
    checkInd = 1;
    while checkInd <= length(inRow),
        if inRow(checkInd),
            nThisRun = 1;
            while (checkInd+nThisRun) <= length(inRow) && inRow(checkInd+nThisRun) == 1,
                nThisRun = nThisRun+1;
            end
            nConsecRow(checkInd:(checkInd+(nThisRun-1))) = nThisRun;
            starts(checkInd:(checkInd+(nThisRun-1))) = checkInd;
            checkInd = checkInd+nThisRun;
        end
        checkInd = checkInd+1;
    end
    out(ir,:) = nConsecRow;
end