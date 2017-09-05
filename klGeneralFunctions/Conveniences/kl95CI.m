function out = kl95CI(in,varargin)

% Set defaults
tail = 'two';
dim = 1;
ci = 95;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t'},
            tail = varargin{varStrInd(iv)+1};
        case {'-c','-ci'},
            ci = varargin{varStrInd(iv)+1};
    end
end

% Transpose if necessary
if dim == 2,
    in = in';
end

% Make sure we're not living in percent land
if ci > 1,
    ci = ci/100;
end

% Sort input
sortIn = sort(in,1);

% Next grab them based on which tails we want
out = nan(2,size(in,2));
switch tail
    case {'two','both'},
        out(1,:) = sortIn(ceil(((1-ci)/2).*size(sortIn,1)),:);
        out(2,:) = sortIn(size(sortIn,1)-floor(((1-ci)/2).*size(sortIn,1)),:);
end
