function [outMat, alignedTimes] = klTroughAlign(inMat,times,tCrit,varargin)

wind = 10;
tInds = find(times == tCrit);
tInc  = nanmean(unique(diff(times)));
minVal = .2;
visualize = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'wind','-w'},
            wind = varargin{varStrInd(iv)+1};
        case {'-v','visualize'},
            visualize = varargin{varStrInd(iv)+1};
    end
end

% Loop through rows of inMat
for ir = 1:size(inMat,1),
    myWind = 1;
    nTries = 0;
    while min(myWind) > minVal/wind && nTries < 100 && (tInds+(wind*nTries)) < (size(inMat,2)-1),
        % Get derivative of the row
        div = diff(inMat(ir,:));
        myWind = abs(div((tInds-wind):(tInds+(wind*nTries+1))));

        % Find the minimum of myWind
        windT = find(myWind == min(myWind),1);
        if isempty(windT), windT = 0; end
            
        if visualize,
            figure();
            plot(inMat(ir,:)'); vline(tInds-wind); vline(tInds+(wind*nTries+1));
            keyboard
        end
        % Get index and save the row - this will be an input to klAlign
        myT(ir) = tInds-wind+windT;
        tempCell{ir,1} = inMat(ir,:);
        nTries = nTries + 1;
    end
    if ~exist('windT','var'),
        windT = find(myWind == min(myWind),1);
        myT(ir) = tInds-wind+windT;
        tempCell{ir,1} = inMat(ir,:);
    end
end

[outMat, tZero] = klAlignv2(tempCell,myT);
tMin = (tZero-1).*(-tInc); tMax = (size(outMat,2)-tZero).*tInc;
alignedTimes = tMin:tInc:tMax;    
    