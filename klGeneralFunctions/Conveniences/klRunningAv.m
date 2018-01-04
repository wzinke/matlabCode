function [out,std] = klRunningAv(in,width)

pad = [nan(1,width),in,nan(1,width)];

runMat = nan(width,length(pad));
for i = 1:width
    runMat(i,:) = [pad(i:end),nan(1,(i-1))];
end
out = nanmean(runMat(:,(width+1):(end-width)),1);
std = nanstd(runMat(:,(width+1):(end-width)),[],1)./sqrt(size(runMat,1));
