function [mnVal,outVals,varVal,outVarVal] = klAssessClustQual(vis,mov,idx)

% Define the function to produce
varRange = @(x) nanmean(range(x,1)./range(nanmean(x,1)));
varRange2 = @(x) nanmean(var(x,1)./var(nanmean(x,1)));
uClusts = unique(idx(~isnan(idx)));
outVals = nan(1,length(uClusts));
outVarVal = nan(1,length(uClusts));
% Loop through
for ic = 1:length(uClusts),
    outVals(ic) = (varRange(vis(idx==uClusts(ic),:)) + varRange(mov(idx==uClusts(ic),:)))/2;
    outVarVal(ic) = (varRange2(vis(idx==uClusts(ic),:)) + varRange2(mov(idx==uClusts(ic),:)))/2; 
end
mnVal = nanmean(outVals);
varVal = nanmean(outVarVal);
