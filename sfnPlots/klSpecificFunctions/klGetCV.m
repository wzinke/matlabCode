%% 15-07-21 function gets coefficient of variation from an input of spiketimes
%
% Can return standard coefficient of variation (default) or local coefficient of variation if the '-type' flag is used
% Cv2 calculated according to equation 4 of:
%      Holt GR, Softky WR, Koch C, Douglas RJ (1996) Comparison of
%      discharge variability in vitro and in vivo in cat visual cortex
%      neurons. J Neurophysiology 75:4806-1814
%
%      

function cv = klGetCV(spiketimes,varargin)

report = 0;
type   = 'standard';

% Decode varargin
varStrInd = cellfun(@ischar,varargin);
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case '-type'
            type = varargin{varStrInd(iv)+1};
    end
end

nSpks = sum(isfinite(spiketimes),2);
% Stop function, return nan if there are fewer than 4 trials of no spikes
if sum(nSpks > 0) < 4,
    if report
        fprintf('Too few spikes to get Cv\n');
    end
    cv = nan;
    return
end

switch type
    case 'standard'
        % Get ISI
        isiMat = diff(spiketimes,1,2);
        % CV = sd/mean
        cv = nanstd(isiMat(:))/nanmean(isiMat(:));
    case 'local'
        % Get ISIs, convert to one long vector
        isiMat = diff(spiketimes,1,2);
        isiVect = nan(1,numel(isiMat));
        for it = 1:size(isiMat,1),
            isiVect = cat(2,isiVect,isiMat(it,:));
        end
        % Initialize some vectors
        cv2 = nan(1,numel(isiVect)-1); 
        % Loop through ISI vector - 
        for ii = 1:length(isiVect)-1,
            cv2(ii) = (2*abs(isiVect(ii+1)-isiVect(ii)))/(isiVect(ii+1)+isiVect(ii));
        end
        cv = nansum(cv2)./(sum(~isnan(cv2))-1);
end
