function [outMat, alignedTimes] = klTroughAlignv5(inMat,times,tCrit,varargin)

wind = 10;
tInds = find(times == tCrit);
tInc  = nanmean(unique(diff(times)));
minVal = .15;
visualize = 0;
tMax = 100;
tMin = -30;
tMax2 = 150;
tMin2 = -80;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'wind','-w'},
            wind = varargin{varStrInd(iv)+1};
        case {'-v','vis'},
            visualize = varargin{varStrInd(iv)+1};
    end
end

% Loop through rows of inMat
for ir = 1:size(inMat,1),
    
    minInd = find(inMat(ir,:) == min(inMat(ir,:)));
    maxInd = find(inMat(ir,:) == max(inMat(ir,:)));
    
    ts = [minInd,maxInd];
    absInd = abs(ts) == min(abs(ts));
    
    if ~isempty(absInd)
        myT(ir) = ts(absInd);
    else
        myT(ir) = tInds;
    end
    
    tempCell{ir,1} = inMat(ir,:);
    
    if visualize,
        figure()
        plot(times,inMat(ir,:));
        vline(times(myT(ir)));
        keyboard
    end
end

[outMat, tZero] = klAlignv2(tempCell,myT);
tMin = (tZero-1).*(-tInc); tMax = (size(outMat,2)-tZero).*tInc;
alignedTimes = tMin:tInc:tMax;    
    