function fnScore = klGetFNv1(spkClust,nzClust,varargin)

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

% Get k nearest neighbors for each noise event
k = ceil(nnPerc.*size(spkClust,1));
isNFA = nan(size(nzClust,1),1);
for inz = 1:size(nzClust,1),
    simNz = [simMat(1:(inz+size(spkClust,1)),(inz+size(spkClust,1)))',simMat((inz+size(spkClust,1)),(inz+size(spkClust,1)+1):end)];
    [sortNz,sortInd] = sort(simNz,'descend');
    knnInds = sortInd(1:k);
    if any(isnan(simNz(knnInds))),
        tmpInds = nan(1,k);
        tmpInds(1:sum(~isnan(simNz(knnInds)))) = knnInds((sum(isnan(simNz(knnInds)))+1):end);
        tmpInds((sum(~isnan(simNz(knnInds)))+1):k) = sortInd((k+1):(k+sum(isnan(simNz(knnInds)))));
        knnInds = tmpInds;
    end
    nSimSpks(inz) = sum(ismember(knnInds,[1:size(spkClust,1)]));
    isNFA(inz) = nSimSpks(inz)/k >= faCut;
end
nfa = sum(isNFA);
fnScore = nfa/(nfa+size(spkClust,1));


%% Similarity function time
function thisSim = getSim(dist,lamb,dNorm)
%     dist = @(x,y) nansum((x-y).^2);
    thisSim = exp(-(dist.*lamb./dNorm));
end

end