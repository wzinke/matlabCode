%% klGetType takes as input spike times and a Task structure with fields for
%  events in the task to test against. Outputs various comparison measures.
%
%  Inputs:
%     - spks: mxn matrix where m = number of trials and n = maximum number
%       of spikes within one trial. Columns contain spike times during the
%       appropriate trial, relative to the beginning of the trial (at
%       least, trials are aligned consistently). Extra entries are nans
%       (e.g., spks = [2 4 6 7 nan; 3 4 5 8 10; 1 3 nan nan nan];)
%     - Task: structure containing fields for events within the trial.
%       Currently, expected fields are Correct, ErrorTone, StimOnset, SRT,
%       Reward Tone, and Reward. These times are expected to be aligned to
%       the same time as spks.
%
%  Outputs:
%     - pOut: Uncorrected p values for the various comparisons made within
%       this function.
%     - compTypes: Names of the tests performed
%     - type: "vis", "mov", or "vismov" - cell classification as determined
%       by uncorrected p values (may not be valid?)
%     - vmi: Visuomotor Index as defined by Cohen et al 2009
%     - spkRate: mxn matrix where m = number of trials  and n = number of
%       conditions during which spike rates were measured.
%     - rateCols: Names corresponding to the intervals of the columns of
%       spkRate (i.e., length(rateCols) == size(spkRate,2)
%
%  Created by Kaleb Lowe June 2015
%  Last edited by Kaleb Lowe 8/15/15

function [type] = klGetTypeAlt(spks,Task,varargin)

% Set defaults
minDur = 10;
meanType = {'mean','median'};
sdThresh = 3;
minSpks = 3.*size(spks,1);

% Decode varargin
if nargin > 2
    % Check that there are enough input arguments, if optional flags are
    % set. None of the following are really used.
    if length(varargin{1}) ~= length(varargin), error('Not enough input arguments'); end
    for iv = 2:length(varargin{1})
        switch varargin{1}(iv)
            case 'c'
                EV          = varargin{iv};
            case 't'    
                thresh      = varargin{iv};
            case 'p'
                doPlot      = varargin{iv};
            case 'm'
                meanType    = varargin{iv};
        end
    end
end

if ~exist('EV','var'),
    load cosmanCodes;
end

if sum(isfinite(spks(:))) < minSpks,
    type = {'none','none'};
    return;
end

%% Define intervals for tests
% These windows are defined relative to the appropriate event. For example,
% postSaccWind = [0, 100] tests the first 100ms after the saccade (i.e., the
% time specified by Task.SRT)

preTrialWind    = [-200, -100];
preVisWind      = [-100, 0];
visTransWind    = [50, 150];
visSustWind     = [150, 250];
prePreSaccWind  = [-200, -100];
preSaccWind     = [-100, 0];
postSaccWind    = [0, 100];
postToneWind    = [0, 100];
preErrWind      = [-100, 0];
postErrWind     = [0, 100];
postRewdWind    = [0, 100];

% Cut post-sacc spikes from vSpks
vSpks = spks;
for ii = 1:size(vSpks,1),
    vSpks(ii,vSpks(ii,:) > Task.SaccEnd(ii)) = nan;
end

%% Get SDF aligned on stimulus and srt
[vSDF, vSDFTimes] = klSpkRatev2(vSpks,'-q',1);
[mSDF, mSDFTimes] = klSpkRatev2(spks-repmat(Task.GoCue + Task.SRT,1,size(spks,2)),'-q',1);
% [mSDF, mSDFTimes] = klSpkRatev2(spks-repmat(Task.SaccEnd,1,size(spks,2)),'-q',1);

%% Get baseline mean and sd for Z score
blMean = nanmean(nanmean(vSDF(Task.Correct == 1,vSDFTimes < preTrialWind(2) & vSDFTimes >= preTrialWind(1)),1),2);
blStd = nanstd(nanmean(vSDF(Task.Correct == 1,vSDFTimes < preTrialWind(2) & vSDFTimes >= preTrialWind(1)),1),[],2);
blMedian = nanmedian(nanmean(vSDF(Task.Correct == 1,vSDFTimes < preTrialWind(2) & vSDFTimes >= preTrialWind(1)),1),2);
blMad = mad(nanmean(vSDF(Task.Correct == 1,vSDFTimes < preTrialWind(2) & vSDFTimes >= preTrialWind(1)),1),1,2);

%% Get sustained visual mean and SD
vsMean = nanmean(nanmean(vSDF(Task.Correct == 1,vSDFTimes >= visSustWind(1) & vSDFTimes <= visSustWind(2)),1),2);
vsSTD = nanstd(nanmean(vSDF(Task.Correct == 1,vSDFTimes >= visSustWind(1) & vSDFTimes <= visSustWind(2)),1),[],2);

for it = 1:length(meanType)
    switch meanType{it}
        case 'mean'
            vz = normData(vSDF,'-z',[blMean,blStd]);
            mz = normData(mSDF,'-z',[blMean,blStd]);
            mzv = normData(mSDF,'-z',[vsMean,vsSTD]);
        case 'median'
            vz = normData(vSDF,'-z',[blMedian,blMad]);
            mz = normData(mSDF,'-z',[blMedian,blMad]);
    end

    % Make threshold sdThresh*SD... Latencies in papers have been defined as the times
    % at which the SDF exceeds 2*SD  *but continues to exceed 6* 
    nConsecVT = klGetConsecutive(abs(nanmean(vz(:,vSDFTimes > visTransWind(1) & vSDFTimes < visTransWind(2)),1)) > sdThresh);
    nConsecVS = klGetConsecutive(abs(nanmean(vz(:,vSDFTimes > visSustWind(1) & vSDFTimes < visSustWind(2)),1)) > sdThresh);
    nConsecPS = klGetConsecutive(abs(nanmean(mz(:,mSDFTimes > preSaccWind(1) & mSDFTimes < preSaccWind(2)),1)) > sdThresh & ...
        abs(nanmean(mzv(:,mSDFTimes > preSaccWind(1) & mSDFTimes < preSaccWind(2)),1)) > sdThresh);

    isVis = 0; isMov = 0;
    if any(nConsecVT > minDur) || any(nConsecVS > minDur),
        isVis = 1;
    end
    if any(nConsecPS > minDur)
        isMov = 1;
    end

    if isVis && isMov,
        type{it} = 'vismov';
    elseif isVis && ~isMov,
        type{it} = 'vis';
    elseif ~isVis && isMov,
        type{it} = 'mov';
    else
        type{it} = 'none';
    end
    % keyboard
end
end
