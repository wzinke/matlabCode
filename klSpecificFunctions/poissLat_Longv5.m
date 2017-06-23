%% This function takes a matrix of spike times (padded with nans) and returns the poisson-based latency
%
%  Good (seemingly) working version named "poissLat_Longv5.m"...
%     Copied to klGeneralFunctions under the name "klPoissLatv1.m" on
%     8/3/15
%
%  This file aims to be a (hopefully overly) comprehensive, step-by-step
%  description of the poisson spike train analysis from Hanes et al 1995.
%  Perhaps these steps will be compressed in a later version, but here, no
%  shortcuts are expected to be taken. This laborious approach aims to help
%  me (and future users as well) understand the process in this method of
%  latency determination and rationale.
%
%
%  "Poisson spike train analysis determines how improbable it is that the
%  number of action potentials within a specific time interval is a chance
%  occurrence."

function [modLat, activeTimes] = poissLat_Longv5(spiketimes,varargin)
% Set some defaults
sigAlph = .005;
checkAlph = .005;
visualize = 0;
stimTime  = 0;
fixTime   = -500;
minChkSpks = 3;
% nCheckSpks = [min([3,length(;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case '-rwd'
            rwdTimes = varargin{varStrInd(iv)+1};
        case '-vis'
            visualize = varargin{varStrInd(iv)+1};
        case '-stim'
            stimTime = varargin{varStrInd(iv)+1};
            fixTime  = stimTime - 500;
        case '-minSpks'
            minChkSpks = varargin{varStrInd(iv)+1};
    end
end

% "The determination of the mean discharge rate, r, is crucial...Thus, r
% was defined as the mean discharge rate of the entire spike train being
% analyzed...we analyzed the spike train recorded during each trial..."
% "The mean discharge rate, r, during the trial being analyzed was
% determined as the number of spikes in the interval from the fixation of
% the central spot to the time the reward was given divided by the duration
% of that interval."
% 
% WZ's data extraction/formatting aligns all spikes to stim onset and
% includes 500ms of baseline (i.e., the fixation window). Therefore our
% to find the spikes between fixation and reward, all spikes are after
% fixation. However, the reward time isn't as straightforward.
% 
% 7-29-15 edit: Turns out post-fixation isn't as straightforward as I
% originally thought. That is, if the spiketimes are realigned on a
% different event (e.g., SRT) then not all fixation times will be -500.
% This being the case, I made a variable "fixTime" that carries this
% information. It defaults to -500 and is set with the other defaults
% above, but it needs to be a column vector with the same number of rows as
% "spiketimes". Same again for "stimTime"
if length(fixTime) == 1,
    fixTime = repmat(fixTime,size(spiketimes,1),1);
elseif length(fixTime) ~= size(spiketimes,1),
    error('Size mismatch: length(fixTime) ~= size(spiketimes,1)');
end
if length(stimTime) == 1,
    stimTime = repmat(stimTime,size(spiketimes,1),1);
elseif length(stimTime) ~= size(spiketimes,1),
    error('Size mismatch: length(stimTime) ~= size(spiketimes,1)');
end

% This gives two options if we're trying to maximize modularity of this
% script. Those options are:
%     1. Giving the vector of reward times as an optional input argument
%     2. Assuming the last spike time in that trial is the reward time
%
% The former option is more appropriate for full analysis replication, the
% second option allows everything to be calculated with one input. Both
% options are handled in the following section:

if ~exist('rwdTimes','var'),
    rwdTimes = nan(size(spiketimes,1),1);
    for it = 1:size(spiketimes,1),
        lastSpkTime  = find(~isnan(spiketimes(it,:)),1,'last');
        if ~isempty(lastSpkTime)
            rwdTimes(it) = spiketimes(it,find(~isnan(spiketimes(it,:)),1,'last'));
        end
    end
end


% Initialize output variables
trMeanRate  = nan(size(spiketimes,1),1);
activeInds  = nan(size(spiketimes,1),2);
activeTimes = nan(size(spiketimes,1),2);
modLat      = nan;

% The next steps will need to be done for all trials, so let's start a
% trial loop
for it = 1:size(spiketimes,1),
    % "The program then indexed to the first spike after target presentation
    % and advanced through the spike train until finding the first two
    % consecutive spikes that had a mean discharge rate >= to r.
    
    % Get r
    thisTrSpks     = spiketimes(it,~isnan(spiketimes(it,:)));
    % If there are no non-nan spikes this trial, continue to the next trial
    if isempty(thisTrSpks), continue; end
    
    trMeanRate(it) = length(thisTrSpks)/(rwdTimes(it)-fixTime(it));
%     trMeanRate(it) = sum(thisTrSpks < stimTime(it))/(stimTime(it)-fixTime(it));

    % Index to first spike after target presentation
    firstSpkInd = find(thisTrSpks >= stimTime(it),1);
    % Because the algorithm requires three spikes (firstSpkInd and the next
    % two) for initial SI calculation, make sure that firstSpkInd + 2 <
    % length(thisTrSpks)
    if isempty(firstSpkInd), continue; end
    if (firstSpkInd+2) > length(thisTrSpks), continue; end    
    
    % Advance through spike train until two consecutive spikes have mean
    % discharge rate >= r (i.e., 2 spks / (t2-t1) >= trMeanRate(it))
    % Initialize variable for while loop: thisPairRate
    startPair    = firstSpkInd;
    thisPairRate = 2/(thisTrSpks(firstSpkInd+1)-thisTrSpks(firstSpkInd));
    
    % Our while loop advances through the train (incrementing startPair)
    % and checking whether the distance between "startPair" and
    % "startPair+1" gives a rate greater than the trial mean. It also exits
    % the while loop when "startPair+1" is a valid incex for thisTrSpks
    while thisPairRate < trMeanRate(it) && (startPair+1) < length(thisTrSpks),
        startPair = startPair + 1;
        thisPairRate = 2/(thisTrSpks(startPair+1)-thisTrSpks(startPair));
    end
    
    % Because the while loop will exit when the whole train has been
    % examined, regardless of whether a putative burst has been found or
    % not, take a moment here to confirm which of those options caused the
    % exit
    if thisPairRate < trMeanRate(it),
        continue;
    end
    
    % "The time between these two spikes was the inital value of T. Then,
    % the next spike was indexed and the ISI between this spike and the
    % previous spike was added to T. The corresponding SI was then
    % calculated. Successive spikes were then indexed and their ISIs were
    % added to T until the end of the spike train. The SI was calculated
    % after each addition... The spike at the end of hte interval T with
    % the maximum SI value was defined as the end of the putative burst"
    
    % That is, let T be the time between thisTrSpks(startPair) and
    % thisTrSpks(startPair + 2). Looking ahead, because we will be adding
    % spikes incrementally, we will let that 2 be a variable,
    % "burstLenAdd", with a value of 1 (startPair + 1 + burstLenAdd)
    burstLenAdd = 1;
    
    % Advance through the train with a while loop...
    forSI = [];
    while (firstSpkInd+1+burstLenAdd) <= length(thisTrSpks),
        forSI = cat(1,forSI,getSI(2+burstLenAdd,(thisTrSpks(firstSpkInd+1+burstLenAdd)-thisTrSpks(firstSpkInd))*trMeanRate(it)));
        burstLenAdd = burstLenAdd+1;
    end
    % Now find the index (corresponding to burstLenAdd value) of the
    % maximized surprise index
    maxLenAdd = find(forSI == max(forSI));
    
    % At least during development, go to debug mode if more than one value
    % is the maximum SI
    if length(maxLenAdd) > 1,
        %keyboard
        % Don't know what the right approach here is, but let's say we go
        % with the latest of the options to make sure we catch the whole
        % burst?
        maxLenAdd = max(maxLenAdd);
    end
    burstEndInd = firstSpkInd+1+maxLenAdd;
    
    % "Next, the program indexed to the last spike in the spike train, and
    % the SI was calculated for the time interval T from the last spike to
    % the first spike after target presentation. The program then removed
    % spikes from the beginning of the spike train until reaching the spike
    % defining the end of the burst [burstEndInd]. The spike at which the
    % SI was maximized was defined as the beginning of the putative burst.
    % If the SI from the beginning of hte burst to the end of the burst was
    % not significant (p < .005 [?]) the trial was defined as having no
    % burst."
    
    % That is, let's loop from startSpk to burstEndInd and calculate the SI
    % from our loop counter to the last spike in the train (not just
    % burstEndInd). This will give us "revSpkInd" which is the index of the
    % first spike in the burst that maximizes SI.
    
    revSI = [];
    for is = firstSpkInd:burstEndInd,
        revSI = cat(1,revSI,getSI(length(thisTrSpks)-(is-1),(thisTrSpks(end)-thisTrSpks(is))*trMeanRate(it)));
    end
    maxRevSI  = find(revSI == max(revSI));
    if length(maxRevSI) > 1,
        %keyboard
        % Don't know what the right approach here is, but let's say we go
        % with the earliest of the options to make sure we catch the whole
        % burst? (Mirroring the strategy from earlier)
        maxRevSI = min(maxRevSI);
    end
    burstStartInd = maxRevSI+firstSpkInd-1;
    
    % Now to check that this burst is significant...
    % If not, no point in checking for the beginning/end of activity, so
    % let's just continue to the next trial. Let's also continue if the
    % either: the "start of the burst" is after or equal to the "end of the burst"
    
    if isempty(burstStartInd) || isempty(burstEndInd), continue; end
    if burstStartInd >= burstEndInd, continue; end
    [checkSI, checkP] = getSI(burstEndInd-(burstStartInd-1),trMeanRate(it)*(thisTrSpks(burstEndInd)-thisTrSpks(burstStartInd)));
    if isempty(checkP) || checkP >= sigAlph, continue; end
    
    % Now we have indices for the beginning and end of bursts (and bursts
    % only! - for example a test case listed this burst as being 2ms long
    % because of two very tightly grouped spikes. However, this doesn't
    % mean that activation doesn't spread beyond these two consecutive
    % spikes. 
    
    % "To determine the beginning of activation, the SI first was
    % calculated from the end of the burst to the beginning of the burst
    % ["checkSI"] and then spikes before the beginning of the burst were
    % added in succession until the SI fell below the desired significance
    % level
    
    % To deal with the leniency issue mentioned below, try the following
    % threshold method.
    nCheckSpks = max([minChkSpks,ceil(length(thisTrSpks)*.1)]);
    
    checkStartInd = burstStartInd;
    for is = checkStartInd:-1:1,
        % This criterion seems to be too generous. I want to re-check using
        % the two latest spikes plus the one to be added
%         [~,checkP] = getSI(burstEndInd-(is-1),trMeanRate(it)*(thisTrSpks(burstEndInd)-thisTrSpks(is)));
        [~,checkP] = getSI(nCheckSpks,trMeanRate(it)*(thisTrSpks(min([is+(nCheckSpks+1),length(thisTrSpks)]))-thisTrSpks(is)));
        if checkP > checkAlph, break; end;
    end
    % Again, since this is a limited loop it would hypothetically be
    % possible for the activation to go all the way back to the first spike
    % in the trial. In any case, checkP holds the most recently calculated
    % p value for spikes added to the beginning of the burst. So if checkP
    % is still < checkAlph, then activation starts at index of 1.
    % Otherwise, the activation stats at index of the most recent is + 1
    % (is+1 because is failed to be significant).
    if checkP <= checkAlph, postCheckStartInd = 1; else postCheckStartInd = is+1; end;
    
    % Now let's do the reverse!
    % "To determine the end of activation, the SI was first calculated from
    % the end of the burst to the beginning of the burst [still checkSI],
    % and then spikes after the end of the burst were added in succession
    % until SI fell below the desired significance level
    checkEndInd = burstEndInd;
    for is = burstEndInd:length(thisTrSpks),
        % This criterion seems to be too generous. I want to re-check using
        % the two latest spikes plus the one to be added
%         [~,checkP] = getSI(is-(burstStartInd-1),trMeanRate(it)*(thisTrSpks(is)-thisTrSpks(burstStartInd)));
        [~,checkP] = getSI(nCheckSpks,trMeanRate(it)*(thisTrSpks(is)-thisTrSpks(max([is-(nCheckSpks-1),1]))));
        if checkP > checkAlph, break; end;
    end
    % As before, the burst could go through the end of the spike train. If
    % this is the case, checkP will still be < checkAlph and the end of
    % activation will be the length(thisTrSpks). Otherwise, end of
    % activation = is-1 (because is was the first rejected spike).
    if checkP <= checkAlph, postCheckEndInd = length(thisTrSpks); else postCheckEndInd = is-1; end
    
    activeInds(it,:)  = [postCheckStartInd,postCheckEndInd];
    activeTimes(it,:) = thisTrSpks(activeInds(it,:));
    
    turnOffFlag = 0;
%     if activeTimes(it,1) < 0 && ~visualize,
%         visualize =1;
%         turnOffFlag = 1;
%     end
%     
    if visualize,
        figure()
        scatter(thisTrSpks,repmat(1,1,length(thisTrSpks)),'r*');
        vline(activeTimes(it,1)); vline(activeTimes(it,2));
        %keyboard
    end
    
    if turnOffFlag,
        visualize = 0;
        turnOffFlag = 0;
    end
    
end
modLat = poissMode(activeTimes(:,1));
    
    function [SI, prob] = getSI(n,rT,type)
        % Type allows us to set the -log(P) logarithm to be natural log or
        % log10 (log10 by default).
        if nargin < 3,
            type = 'base10';
        end

        % Let's define the expression for poisson probability. As per Hanes et al
        % 1995:
        %   P = e^(-rT)*EXP where EXP is the sum from n to infinity of ((rT)^i)/i!
        %   Since this is the probability of  N OR MORE spikes, this should be
        %   identical to calculating 1-poisscdf((n-1),rT), a much easier value to
        %   code in MATLAB. 

        % as per the above:
        prob = 1-poisscdf((n-1),rT);

        % Now because SI = -log(P),
        switch type
            case 'base10'
                SI = -(log10(prob));
            case 'natural'
                SI = -(log(prob));
            otherwise
                % Default to base 10
                SI = -(log10(prob));
        end
    end

    function modLat = poissMode(startTimes)
        %% Start mode analysis as described in Thompson et al 1996
        % "First, the times of the beginnings of responses across all trials were
        % sorted in ascending order.
%         activeStarts = activeTimes(~isnan(activeTimes(:,1)),1);
        % Cut out nans
        startTimes(isnan(startTimes)) = [];
        [sortStarts, sortInds] = sort(startTimes,'ascend');

        % "Next, an integer J was selected as the window size over which the
        % distribution function was smoothed. Reliable estimates of the mode
        % resulted when J was set as 1/4 of the total number of periods of
        % significant activations (N) but no less than 3. Next, for every i =
        % 1,...,N-J the mode test statistic (p) was estimated by the equation:
        %  p[((t(i) + t(i+J))/2] = J/N(t(i+J)-ti)

        J = ceil(length(sortStarts)/4); 
        if J < 3
            J = 3; % Set to 3, but this may be a bad sign for mode estimation...
        end
        pJ = [];
        for i = 1:(length(sortStarts)-J),
        %     pJ((sortStarts(i)+sortStarts(i+J))/2) = J/(length(sortStarts)*(sortStarts(i+J)-sortStarts(i)));
            pJ = cat(2,pJ,J/(length(sortStarts)*(sortStarts(i+J)-sortStarts(i))));
        end

        % The time of (t(i)+t(i+J))/2 that generates the largest value of p was the
        % estimated mode and represented the visual latency of the cell
        maxInd = find(pJ == max(pJ));
        modLat = (sortStarts(maxInd(1))+sortStarts(maxInd(1)+J))/2;
    end
end
    
    
    


