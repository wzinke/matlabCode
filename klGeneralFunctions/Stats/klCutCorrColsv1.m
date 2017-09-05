function [cutMat, keepCols] = klCutCorrColsv1(inMat,varargin)

% Set defaults
startCol = 1;
pCrit    = .1;
rCrit    = .5;
type     = 'r';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-s'},
            startCol = varargin{varStrInd(iv)+1};
        case {'-p'},
            pCrit   = varargin{varStrInd(iv)+1};
        case {'-r'},
            rCrit = varargin{varStrInd(iv)+1};
        case {'-t'},
            type = varargin{varStrInd(iv)+1};
    end
end

keepCorr = [];
cutMat = inMat;
keepCols = 1:size(inMat,2);
for iFact1 = startCol:(size(cutMat,2)-1),
    sigCorrCut = [];
    for iFact2 = (iFact1+1):(size(cutMat,2)),
        clear fact1 fact2
        fact1 = inMat(:,iFact1);
        fact2 = inMat(:,iFact2);
        [r,p] = corr(fact1(~isnan(fact1) & ~isnan(fact2)),fact2(~isnan(fact1) & ~isnan(fact2)));
%         [r,p] = corr(fact1,fact2);
        if isnan(r),
            keyboard
        end
        switch type
            case {'p'},
                if p < pCrit
                    sigCorrCut(length(sigCorrCut)+1) = iFact2;
                    keepCols = [keepCols(1:(iFact2-1)),keepCols((iFact2+1):end)];
        %             keepCorr(length(keepCorr)+1) = startCol+(iFact1-1);
    %                 cutMat = cutMat(:,~ismember(1:size(cutMat,2),iFact2));
                end
            case {'r'},
                if abs(r) > rCrit
                    sigCorrCut(length(sigCorrCut)+1) = iFact2;
%                         keepCols = [keepCols(1:(iFact2-1)),keepCols((iFact2+1):end)];
        %             keepCorr(length(keepCorr)+1) = startCol+(iFact1-1);
    %                 cutMat = cutMat(:,~ismember(1:size(cutMat,2),iFact2));
                end
        end

    end
    keepCols = keepCols(~ismember(1:size(cutMat,2),sigCorrCut));
    cutMat = inMat(:,keepCols);
end
