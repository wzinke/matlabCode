function [outK,outIDs,outGap,outErr,scores,coeffs,eigVect] = klPlotSortCombos(wvs,dimRedType,sortType,figNum,vis,varargin)

colors = 'rgbcmky';

varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'ndims'},
            nDims = varargin{varStrInd(iv)+1};
    end
end

switch dimRedType
    case 'pca',
        [coeffs,scores] = pca(wvs,'Centered','off');
        eigVect = nan;
    case 'lpp',
        [scores, eigVect] = klLPPv1(wvs);
        coeffs = nan;
end


if ~exist('nDims','var'),
    nDims = size(scores,2);
end

switch sortType,
    case 'kmeans',
        [outK, outIDs, outGap, outErr] = klAutoSortv3a(scores(:,1:nDims),'-p',1,'-t','pca','perc',1);
    case 'lsc'
        [outK, outIDs, outGap, outErr] = klLSCwGap(scores(:,1:nDims),'-p',min([ceil(size(wvs,1)/2),1000]),'perc',1);
end

if vis,
    if nargin < 4,
        figure();
    else
        figure(figNum);
    end

    if outK == 1,
        subplot(2,1,1);
    else
        subplot(outK,outK,1:(outK*(outK-1)));
    end

    for iw = 1:outK,
        scatter(scores(outIDs(:,outK)==iw,1),scores(outIDs(:,outK)==iw,2),[],colors(iw));
        hold on;    
    end

    % figure();
    for iw = 1:outK,
        if outK==1,
            subplot(2,1,2);
        else
            subplot(outK,outK,(outK*outK)-(iw-1));
        end
        plot(wvs(outIDs(:,outK)==iw,:)',colors(iw)); hold on;
    end
end