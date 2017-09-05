function [outMat, alignedTimes, goodShift, alignInd] = klTroughAlignv5(inMat,times,tCrit,varargin)

wind = 10;
minVal = .2;
visualize = 0;
maxShft = 40;

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

% Get the index of the alignment time
tZero = find(times==tCrit);

% Convert to Z to avoid amplitude problems
zMat = (inMat-repmat(nanmean(inMat,2),1,size(inMat,2)))./repmat(nanstd(inMat,[],2),1,size(inMat,2));

% Loop through rows of inMat
goodShift = ones(size(inMat,1),1);
alignInd = ones(size(inMat,1),1).*tZero;
for ir = 1:size(zMat,1),
    
    if sum(isnan(inMat(ir,:)))==size(inMat,2),
        tempCell{ir,1} = nan;
        alignInd(ir) = nan;
        continue;
    end
    
    % Get extrema
    [~,imax,~,imin] = extrema(zMat(ir,:));
    allInds = [imax,imin];
    
    % Get the two closest extrema to "time zero"
    [~, sortIndInd] = sort(abs(allInds-tZero));
    extInd = allInds(sortIndInd(1:2));
    
    % Pick the one with the bigger magnitude
    alignInd(ir) = extInd(abs(zMat(ir,extInd))==max(abs(zMat(ir,extInd))));
    tempCell{ir,1} = inMat(ir,:);
    
    if abs(alignInd(ir)-tZero) > maxShft
        alignInd(ir) = tZero;
        goodShift(ir) = 0;
    end
    
    
        if visualize,
            figure();
            plot(times,inMat(ir,:)'); 
            v(1) = vline(times(extInd(1))); 
            v(2) = vline(times(extInd(2)));
            set(v(abs(zMat(ir,extInd))==max(abs(zMat(ir,extInd)))),'color','r');
            
            keyboard
        end
end

[outMat, tZero] = klAlignv2(tempCell,alignInd);
tInc = nanmean(diff(times));
tMin = (tZero-1).*(-tInc); tMax = (size(outMat,2)-tZero).*tInc;
alignedTimes = tMin:tInc:tMax;    
    