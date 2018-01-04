function out = klRunningAvv2(in,width)

pad = [nan(size(in,1),width),in,nan(size(in,1),width)];

runMat = nan(size(pad,1),size(pad,2),2*width+1);
for i = (-width):width
    if i < 0
        tmp = [nan(size(pad,1),(abs(i))),pad(:,1:(end-abs(i)))];
    elseif i > 0
        tmp = [pad(:,i:end),nan(size(pad,1),(i-1))];
    else
        tmp = pad;
    end
    runMat(:,:,(i+width+1)) = tmp;
end
out = nanmean(runMat(:,(width+1):(end-width),:),3);