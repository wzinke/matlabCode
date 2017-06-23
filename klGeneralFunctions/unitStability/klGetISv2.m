function [isoScore, simMat] = klGetISv2(spkClust,nzClust,varargin)

% Set defaults
silent = 1;
distType = 'euc';
lamb = 10;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'silent','-q'},
            silent = varargin{varStrInd(iv)+1};
        case {'-d'},
            distType = varargin{varStrInd(iv)+1};
        
    end
end


% Concatenate all events
allEvents = [spkClust;nzClust];

if ~exist('simMat','var'),
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
    elseif strcmp(distType,'euc'),
        distMat = EuDist2(allEvents,allEvents,1);
    elseif strcmp(distType,'sqeuc'),
        distMat = EuDist2(allEvents,allEvents,0);
    else
        eucFun = @(x,y) sqrt(nansum((repmat(x,size(y,1),1)-y).^2,2));
        distVect = pdist(allEvents,eucFun);
        distMat = squareform(distVect);
        
%         distMat = nan(size(allEvents,1));
%         for is = 1:size(allEvents,1),
%     %         for iis = (is):size(allEvents,1),
%     %             if iis ~= is,
%     %                 distMat(is,iis) = klGetDist(allEvents(is,:),allEvents(iis,:),'-t',distType);
%     %             end
%     %         end
%             thisMult = allEvents-repmat(allEvents(is,:),size(allEvents,1),1);
%             thisMult(isnan(thisMult)) = 0;
%             distMat(:,is) = diag(thisMult*thisMult');
% 
% 
%         end
        for i = 1:size(distMat,1),
            distMat(i:end,i) = deal(nan);
        end
    end
    if ~silent, fprintf('Distmat created in %s\n',printTiming(distTic)); end

    spkDist = distMat(1:size(spkClust,1),1:size(spkClust,1));
    d0 = nanmean(spkDist(:));
    simMat = getSim(distMat,lamb,d0);
end

if isempty(nzClust),
    isoScore = nan;
    return;
end

if ~silent, fprintf('Getting all P(X)...'); pxTic = tic; end
% Now get all P(X) where X is an event in the spike cluster
PX = zeros(size(spkClust,1),1);
for ix = 1:size(spkClust,1),
    % Get sum distance to other clusters
% 	sumSim = nansum(simMat(ix,:));
    sumSim = nansum([simMat(1:ix,ix)',simMat(ix,(ix+1):end)]);
    
    % Now for all Y's, get PxY
%     PX(ix) = nansum(simMat(ix,1:size(spkClust,1))./sumSim);
    PX(ix) = nansum([simMat(1:ix,ix)',simMat(ix,(ix+1):size(spkClust,1))]./sumSim);
end
if ~silent, fprintf('Done in %s\n',printTiming(pxTic)); end

isoScore = nanmean(PX);


%% Similarity function time
function thisSim = getSim(dist,lamb,dNorm)
%     dist = @(x,y) nansum((x-y).^2);
    thisSim = exp(-(dist.*lamb./dNorm));
end

end