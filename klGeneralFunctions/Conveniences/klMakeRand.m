function randVals = klMakeRand(n,params1,params2,varargin)

% Set defaults
dists = ones(1,length(params1));

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)},
        case {'-d','dists','dist'},
            dists = varargin{varStrInd(iv)+1};
    end
end

% Initialize output
randVals = nan(n,length(params1));

% Get gaussian distributions
isGauss = dists==1;
randTmp = randn(n,sum(isGauss));
% Here, params2 = SD and params1 = mu
randTmp = (randTmp.*repmat(params2(isGauss),n,1))+repmat(params1(isGauss),n,1);
randVals(:,isGauss) = randTmp;
clear randTmp;

% Get uniform distributions
isUnif = dists==2;
if any(isUnif),
    randTmp = rand(n,sum(isUnif));
    % Here, params1 = min, params2=max
    randTmp = (randTmp.*repmat(params2(isUnif)-params1(isUnif),n,1))+repmat(params1(isUnif),n,1);
    randVals(:,isUnif) = randTmp;
end