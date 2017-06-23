function [t2s, tStart] = klSigLat(spiketimes,varargin)

% Set defaults
bWind       = 200;
cWind       = 10;
tstWind     = [50,150];
consecTime  = 30;
nStd        = 2;
rev = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'bWind','-b'}
            bWind = varargin{varStrInd(iv)+1};
        case {'-c','cwind','cWind'}
            cWind = varargin{varStrInd(iv)+1};
        case {'-r','rev','dir','dim'}
            rev = varargin{varStrInd(iv)+1};
        case {'baseVals','baseline','base'}
            blVals = varargin{varStrInd(iv)+1};
    end
end

% Initialize output variables
tStart = nan; %t2s defined and dealt with later

% First get SDF
[sdf,sdfTimes] = klSpkRatev2(spiketimes,'-q',1);

% Get baseline values - across trials then across times. Or pull from input
% values
mnSDF = nanmean(sdf,1);
if ~exist('blVals','var')
    blMean = nanmean(mnSDF(sdfTimes >= -bWind & sdfTimes < 0));
    blStd = nanstd(mnSDF(sdfTimes >= -bWind & sdfTimes < 0));
else
    blMean = blVals(1);
    blStd  = blVals(2);
end

% Check if the test window is an enhancement or a suppression
if rev
    tstMean = nanmean(mnSDF(sdfTimes <= -tstWind(1) & sdfTimes >= -tstWind(2)));
    isSupp  = tstMean < blMean;
else
    tstMean = nanmean(mnSDF(sdfTimes >= tstWind(1) & sdfTimes <= tstWind(2)));
    isSupp  = tstMean < blMean;
end

% Find the first time where mnSDF > nStd*blStd for a minimum of consecTime ms
if isSupp
    nConsec = klGetConsecutive(mnSDF < (blMean-nStd*blStd));
else
    nConsec = klGetConsecutive(mnSDF > (blMean+nStd*blStd));
end

if rev
    t2s     = sdfTimes(find(nConsec < consecTime & sdfTimes < 0,1,'last')+1);
else
    t2s     = sdfTimes(find(nConsec >= consecTime & sdfTimes >= 0,1));
end

if isempty(t2s), t2s = nan; end
if ~isnan(t2s),
    
    % Now go backwards to test correlation between mnSDF and time. This reveals
    % the time that the increase starts (as per Sato & Schall 2001)
    binCent = t2s;
    thisSet = mnSDF(sdfTimes >= (binCent-cWind) & sdfTimes <= (binCent + cWind));
    [~,p] = corr(thisSet',(1:length(thisSet))');
    while p < .05 && binCent > (sdfTimes(1)+cWind),
        % binCent > sdfTimes(1)+cWind criterion added to prevent errors
        binCent = binCent-1;
        thisSet = mnSDF(sdfTimes >= (binCent-cWind) & sdfTimes <= (binCent + cWind));
        [~,p] = corr(thisSet',(1:length(thisSet))');
    end
    if p > .05, tStart = binCent+1; else tStart = nan; end
    % Returns tStart = nan if a start time cannot be identified
    % If start time is identified, add 1 to binCent because we want the last
    % significant time bin, not the first non-significant one (working
    % backwards, that is)
end