function [out, outInd] = klClipTails(in,varargin)

% Set defaults
alph = .05;
tails = 'two';

% Decode varargin
varStrInd = find(cellfun('isclass',varargin,'char'));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'tails','tail'}
            tails   = varargin{varStrInd(iv)+1};
        case {'alph','alpha'}
            alph    = varargin{varStrInd(iv)+1};
    end
end

[sortIn, sortInd] = sort(in,'ascend');
lastVal = find(isfinite(sortIn),1,'last');
out = nan(size(in));

switch tails
    case 'upper'
        minInd = 1;
        maxInd = min([length(in),ceil(lastVal*(1-alph))]);
    case 'lower'
        minInd = max([1,floor(lastVal*alph)]);
        maxInd = lastVal;
    case 'two'
        minInd = max([1,floor(lastVal*(alph/2))]);
        maxInd = min([length(in),ceil(lastVal*(1-(alph/2)))]);
    otherwise
        % Default to two-tails
        minInd = max([1,floor(lastVal*(alph/2))]);
        maxInd = min([length(in),ceil(lastVal*(1-(alph/2)))]);
end

out(sortInd(minInd:maxInd)) = in(sortInd(minInd:maxInd));
outInd = sortInd(minInd:maxInd);