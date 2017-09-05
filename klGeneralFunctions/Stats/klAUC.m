%% Estimates area under the curve. Essentially is integrate without having to pass in a function

function auc = klAUC(data,varargin)

% Set defaults
xBins = 1:length(data);
type = 'trap';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case 'x'
            xBins = varargin{varStrInd(iv)+1};
        case {'-t','type'}
            type = varargin{varStrInd(iv)+1};
    end
end

% Make length(data) run in dimension 2
if size(data,2) ~= length(xBins)
    data = data';
end

% Set up difference matrix
x1 = [nan,xBins];
x2 = [xBins,nan];

deltX = x2-x1;
deltX(isnan(deltX)) = 0;

% Get data points by type
switch type
    case {'trap','trapezoid'}
        y = mean([nan,data;data,nan]);
        y(isnan(y)) = 0;
end

% Area = deltX*y
auc = deltX*y';