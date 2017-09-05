function [outMat, alignedTimes] = klTroughAlignv4(inMat,times,tCrit,varargin)

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
    myWind = 1;
    nTries = 0;
    [xmax,imax,xmin,imin] = extrema(inMat(ir,:));
%     imax(xmax < (nanmean(inMat(:,1),1)+nanstd(inMat(ir,:)))) = [];
%     xmax(xmax < (nanmean(inMat(:,1),1)+nanstd(inMat(ir,:)))) = [];
%     imin(xmin > (nanmean(inMat(:,1),1)-nanstd(inMat(ir,:)))) = [];
%     xmin(xmin > (nanmean(inMat(:,1),1)-nanstd(inMat(ir,:)))) = [];
%     
%     imax(xmax < (nanmean(inMat(:,1),1)+nanstd(inMat(:,1)))) = [];
%     xmax(xmax < (nanmean(inMat(:,1),1)+nanstd(inMat(:,1)))) = [];
%     imin(xmin > (nanmean(inMat(:,1),1)-nanstd(inMat(:,1)))) = [];
%     xmin(xmin > (nanmean(inMat(:,1),1)-nanstd(inMat(:,1)))) = [];
    
    diffInds = [imax,imin] - tInds;
    if any(diffInds < 0),
        pre = diffInds(diffInds == max(diffInds(diffInds < 0)));
    end
    if any(diffInds >= 0),
        post = diffInds(diffInds == min(diffInds(diffInds >= 0)));
    end
    if ~exist('pre','var'),
        myT(ir) = post+tInds;
    elseif ~exist('post','var'),
        myT(ir) = pre+tInds;
    else
        slope = inMat(ir,post+tInds)-inMat(ir,pre+tInds);
        isPos = nanmean(inMat(:,tInds)) > 0;

        if (~isPos && slope < 0) || (isPos && slope > 0),
            myT(ir) = post+tInds;
        else
            myT(ir) = pre+tInds;
        end
    end
    
    if ~ismember(myT(ir),[imin,imax]) || myT(ir) > tMax2 || myT(ir) < tMin2,
%         keyboard
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
    