% Function takes as input nxp data matrix "dataMat" where n is the number
% of samples and p is the number of measurements per sample.
% Output is an nxn dissimilarity matrix where n is the number of samples. p

function simMat = klSimilarity(dataMat,varargin)

% Set defaults
type = 'euclidean';
normType = 'range';
normSim  = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'-t','type'}
            type = varargin{varStrInd(iv)+1};
        case {'-n','norm'}
            normSim = varargin{varStrInd(iv)+1};
    end
end

% For starters, let's normalize the columns
switch normType
    case 'range' % Scales 0-1
        minVal = min(dataMat,[],1);
        maxVal = max(dataMat,[],1);
        normMat = (dataMat - repmat(minVal,size(dataMat,1),1))./repmat(maxVal-minVal,size(dataMat,1),1);
    case {'z','std','var'} % Divides by SD
        normMat = dataMat./repmat(nanstd(dataMat,[],1),size(dataMat,1),1);
end

% Loop through the samples and calculate similarity
switch type
    case {'euclidean'} % Euclidean distance (sqrt sum of squared distances)
        simMat = nan(size(dataMat,1));
        % Comparator 1
        for ii = 1:size(dataMat,1),
            % Comparator 2
            for ij = 1:size(dataMat,1),
                % Get the distances
                simMat(ii,ij) = sqrt(nansum((normMat(ii,:)-normMat(ij,:)).^2));
            end
        end
    case {'corr','correlation'} % Correlation measures
        for ii = 1:size(dataMat',1),
            for ij = ii:size(dataMat,1),
                cRho = corr(dataMat(ii,:)',dataMat(ij,:)');
                simMat(ii,ij) = cRho;
                simMat(ij,ii) = cRho;
            end
        end
        
end

if normSim
    maxSim = max(simMat(:));
    minSim = min(simMat(:));
    simMat = (simMat - minSim)./(maxSim-minSim);
end
