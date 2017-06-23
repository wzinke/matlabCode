function [outMat, outX] = klKDFv1(inMat,varargin)

% Set defaults
kernType = 'gauss';
kernWd   = 10;
res      = -1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-k','kern'},
            kernType = varargin{varStrInd(iv)+1};
        case {'-w','width'},
            kernWd = varargin{varStrInd(iv)+1};
        case {'res','-r'},
            res = varargin{varStrInd(iv)+1};
    end
end

% Get resolution
resMat = round(inMat./(10^res));
outX = (10^res):(10^res):max(inMat(:));

% Make matrix of indices
indMat = zeros(size(resMat,1),max(resMat(:)));
for ir = 1:size(indMat,1),
    for ic = 1:sum(~isnan(resMat(ir,:))),
        indMat(ir,resMat(ir,ic)) = indMat(ir,resMat(ir,ic)) + 1;
    end
end

kern = klGetKern('type',kernType,'width',kernWd);
        
for iz = 1:size(inMat,3),
    outMat = conv2(indMat(:,:,iz),kern,'same');
end