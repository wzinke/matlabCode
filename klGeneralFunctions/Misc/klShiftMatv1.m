function [outMat, outT] = klShiftMatv1(inMat,inT,alT)

tInc  = nanmean(unique(diff(inT)));

% Shift time vector
for ir = 1:size(inMat),
    matCell{ir,1} = inMat(ir,:);
    if isnan(alT(ir)), tInds(ir,1) = nan; continue; end
    tInds(ir,1) = find(inT >= alT(ir),1);
    
end

% Align stuff
[outMat,tZero] = klAlignv5(matCell,tInds);
tMin = (tZero-1).*(-tInc); tMax = (size(outMat,2)-tZero).*tInc;
outT = tMin:tInc:tMax;