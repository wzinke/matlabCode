%% 15-07-21 function gets Local Variation from an input of spiketimes
%
% Lv calculated according to equation 2.2 of:
%      Shinomoto S, Shima K, Tanji J (2003) Differences in spiking patterns
%      among cortical neurons. Neural COmput. 15:2823-2842
%
% LvR (revised LV) calculated according to equation 3 of:
%      Shinomoto S, Kim H, Shimokawa T, Matsuno N, Funahashi S, Shima K, Fujita I,
%      Tamura H, Doi T, Kawano K, Inaba N, Fukushima K, Kurkin S, Kurata K,
%      Taira M, Tsutsui K, Komatsu H, Ogawa T, Koida K, Tanji J, et al. (2009)
%      Relating neuronal firing patterns to functional differentiation of cerebral cortex.
%      PLoS Comput Biol 5:e1000433
%      
%      LvR uses a refractory constant of 5ms by default (as per Ardid et al
%      2015)

function lv = klGetLV(spiketimes,varargin)

report = 0;
type   = 'standard';
refract = 5;

% Decode varargin
varStrInd = cellfun(@ischar,varargin);
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case '-type'
            type = varargin{varStrInd(iv)+1};
        case '-refract'
            refract = varargin{varStrInd(iv)+1};
    end
end

nSpks = sum(isfinite(spiketimes),2);
% Stop function, return nan if there are fewer than 4 trials of no spikes
if sum(nSpks > 0) < 4,
    if report
        fprintf('Too few spikes to get Lv\n');
    end
    lv = nan;
    return
end

% Get ISIs, convert to one long vector
isiMat = diff(spiketimes,1,2);
isiVect = nan(1,numel(isiMat));
for it = 1:size(isiMat,1),
    isiVect = cat(2,isiVect,isiMat(it,:));
end
% Initialize some vectors
lvVect = nan(1,numel(isiVect)-1); 
% Loop through ISI vector - 
switch type
    case 'standard'
        for ii = 1:length(isiVect)-1,
            lvVect(ii) = (3*((isiVect(ii)-isiVect(ii+1))^2))/((isiVect(ii)+isiVect(ii+1))^2);
        end
        lv = nansum(lvVect)./(sum(~isnan(lvVect))-1);
    case 'revised'
        for ii = 1:length(isiVect)-1,
            lvVect(ii) = (1-((4*isiVect(ii)*isiVect(ii+1))/((isiVect(ii)+isiVect(ii+1))^2)))*(1+((4*refract)/(isiVect(ii)+isiVect(ii+1))));
        end
        lv = nansum(lvVect).*(3/(length(isiVect)-1));

end
