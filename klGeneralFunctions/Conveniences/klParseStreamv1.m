function [sigMat,tVect,indMat] = klParseStreamv1(sigVect,alignTimes,varargin)

% Set defaults
sigFreq = (1/.983);
% sigFreq = (1000/64.2475);
sigWind = [-500,2500];

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-w'},
            sigWind = varargin{varStrInd(iv)+1};
        case {'-f'},
            sigFreq = varargin{varStrInd(iv)+1};
    end
end

% Get a vector of times
tDiff = (1/sigFreq);
% tDiff = .983;
% tVect = (0:(length(sigVect)-1)).*tDiff;

% Get indices of times
alignInds = round(alignTimes./tDiff);

% Let's make this easier on ourselves and replace the NaNs with 1s for
% now...
alignInds(isnan(alignInds)) = 1;

% Now let's expand this into an index matrix...
startDiff = floor((abs(sigWind(1)/tDiff))*((-1)^(sigWind(1) < 0)));
endDiff = ceil((abs(sigWind(2)/tDiff))*((-1)^(sigWind(2) < 0)));
indMat = bsxfun(@plus,alignInds,startDiff:endDiff);
tVect = (startDiff:endDiff).*tDiff;

% Deal with the old nan situations...
indMat(alignInds==1,:) = 1;
indMat(indMat > length(sigVect)) = 1;

% Get these indices
sigMat = nan(size(indMat));
sigMat(1:numel(sigMat)) = sigVect(indMat(1:numel(indMat)));
sigMat(indMat==1) = nan;

% for i = 1:length(alignInds),
%     sigMat(i,:) = sigVect(indMat(i,:));
% end

% keyboard
