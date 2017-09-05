function fpScore = klGetFPv1(spkClust,nzClust,varargin)

% Set defaults
nnPerc = .05;
faCut = .5;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'sim'},
            simMat = varargin{varStrInd(iv)+1};
        case {'-n','nnperc'}
            nnPerc = varargin{varStrInd(iv)+1};
        
    end
end

if ~exist('simMat','var'),
    % Concatenate all events
    allEvents = [spkClust;nzClust];

    % Get d0 = average Euclidean distance in spike cluster
    if ~silent,
        fprintf('Initializing similarity matrix\n'); distTic = tic;
    end
    if strcmp(distType,'corr'),
        for i = 1:size(allEvents,2),
            goodCols(i) = sum(~isnan(allEvents(:,i))) == size(allEvents,1);
        end

        distMat = 1-corr(allEvents(:,goodCols)',allEvents(:,goodCols)');
        for i = 1:size(distMat,1),
            distMat(i:end,i) = deal(nan);
        end
    else
        distMat = nan(size(allEvents,1));
        for is = 1:size(allEvents,1),
    %         for iis = (is):size(allEvents,1),
    %             if iis ~= is,
    %                 distMat(is,iis) = klGetDist(allEvents(is,:),allEvents(iis,:),'-t',distType);
    %             end
    %         end
            thisMult = allEvents-repmat(allEvents(is,:),size(allEvents,1),1);
            thisMult(isnan(thisMult)) = 0;
            distMat(:,is) = diag(thisMult*thisMult');


        end
        for i = 1:size(distMat,1),
            distMat(i:end,i) = deal(nan);
        end
    end
    if ~silent, fprintf('Distmat created in %s\n',printTiming(distTic)); end

    spkDist = distMat(1:size(spkClust,1),1:size(spkClust,1));
    d0 = nanmean(spkDist(:));
    simMat = getSim(distMat,lamb,d0);
end

% Get k nearest neighbors for each spike event
k = ceil(nnPerc.*size(spkClust,1));
for ispk = 1:size(spkClust,1),
    simSpk = [simMat(1:ispk,ispk)',simMat(ispk,(ispk+1):end)];
    [sortSpk,sortInd] = sort(simSpk,'descend');
    knnInds = sortInd(1:k);
    if any(isnan(simSpk(knnInds))),
        tmpInds = nan(1,k);
        tmpInds(1:sum(~isnan(simSpk(knnInds)))) = knnInds((sum(isnan(simSpk(knnInds)))+1):end);
        tmpInds((sum(~isnan(simSpk(knnInds)))+1):k) = sortInd((k+1):(k+sum(isnan(simSpk(knnInds)))));
        knnInds = tmpInds;
    end
    fpK(ispk,:) = knnInds;
    nSimNz = sum(~ismember(knnInds,[1:size(spkClust,1)]));
    isNFP(ispk) = nSimNz/k >= faCut;
end
nfp = sum(isNFP);
fpScore = nfp/size(spkClust,1);


%% Similarity function time
function thisSim = getSim(dist,lamb,dNorm)
%     dist = @(x,y) nansum((x-y).^2);
    thisSim = exp(-(dist.*lamb./dNorm));
end

end