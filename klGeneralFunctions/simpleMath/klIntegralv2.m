function out = klIntegralv2(in,varargin)

% Set defaults
c = 0;
step = 1;
type = 'trap';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'-c','c','constant'},
            c = varargin{varStrInd(iv)+1};
        case {'-s','step'},
            step = varargin{varStrInd(iv)+1};
        case {'-t','type'},
            type = varargin{varStrInd(iv)+1};
    end
end

out = nan(size(in)); out(1) = c;
switch type
    case 'trap',
        for i = 1:length(in),
            out(i+1) = out(i)+in(
        
out = out.*step;