function [out, outSort] = klCutTails(in,varargin)

% Set defaults
tail = 'two';
alpha = .05;

% Decode varargin
varStrInd = cellfun(@ischar,varargin);
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case 'tail'
            tail = varargin{varStrInd(iv)+1};
    end
end

out = nan(size(in));
for ic = 1:size(in,2)
    [sortIn, sortInds]  = sort(in(:,ic),'ascend');
    noNans              = sortIn(~isnan(sortIn));

    switch tail
        case 'upper',
            cutStart = 1;
            cutEnd   = min([length(noNans),ceil((1-alpha)*length(noNans))]);
        case 'lower',
            cutStart = max([1,floor(alpha*length(noNans))]);
            cutEnd = length(noNans);
        case 'two'
            cutStart = max([1,floor((alpha/2)*length(noNans))]);
            cutEnd   = min([length(noNans),ceil((1-(alpha/2))*length(noNans))]);
    end

    out(sortInds(cutStart:cutEnd),ic)   = in(sortInds(cutStart:cutEnd),ic);
    outSort(cutStart:cutEnd,ic)         = noNans(cutStart:cutEnd);
end    