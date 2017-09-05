function [outN,outC] = klCellHist(inCell)

outC = unique(inCell);
outN = nan(size(outC));
for ic = 1:length(outC),
    outN(ic) = sum(ismember(inCell,outC{ic}));
end