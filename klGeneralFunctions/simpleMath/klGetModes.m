function xOut = klGetModes(x,y,k)

wdLim = 30;

if length(y) ~= length(x),
    y=hist(y,x);
end

[sortY,sortInd] = sort(y,'descend');
sortX = x(sortInd);

tmpY = sortY;
tmpX = sortX;

modeX = sortX(1); tmpX(1) = nan;
modeY = sortY(1); tmpY(1) = nan;
startInd = 2;
nLoops = 0;
while length(modeX) < k && nLoops < length(sortY),
    nextInd = find(~ismember(tmpY,modeY) & ~isnan(tmpY),1);
    if any(abs(modeX-sortX(nextInd)) < wdLim),
        tmpY(nextInd) = nan;
    else
        modeX = cat(2,modeX,sortX(nextInd));
        modeY = cat(2,modeY,sortY(nextInd));
    end
    nLoops = nLoops+1;
end

xOut = modeX;

    