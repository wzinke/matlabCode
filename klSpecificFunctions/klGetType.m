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
%  Last edited by Kaleb Lowe 7/15/15

function [pOut, compTypes, type, vmi, spkRate, rateCols] = klGetType(spks,Task,varargin)

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

%% Define intervals for tests
% These windows are defined relative to the appropriate event. For example,
% postSaccWind = [0, 100] tests the first 100ms after the saccade (i.e., the
% time specified by Task.SRT)

preTrialWind    = [-100, 0];
preVisWind      = [-100, 0];
visTransWind    = [50, 150];
visSustWind     = [150, 300];
prePreSaccWind  = [-200, -100];
preSaccWind     = [-100, 0];
postSaccWind    = [0, 100];
postToneWind    = [0, 100];
preErrWind      = [-100, 0];
postErrWind     = [0, 100];
postRewdWind    = [0, 100];

%% Initialize variables to fill
[preTrial,preVis,visTrans,visSust,prePreSacc,preSacc,postSacc,postTone,postRewd] = deal(nan(1,size(spks,1)));
corrTr = find(Task.Correct);
errTr  = find(~isnan(Task.ErrorTone));

%% Get Rates - find the number of spikes within the appropriate window and divide by that time window (/1000 to make spks/second)
for ii = 1:length(corrTr),
    it = corrTr(ii); 
    preTrial(it)    = sum(spks(it,:) > preTrialWind(1) & spks(it,:) < preTrialWind(2))/(diff(preTrialWind)/1000);
    preVis(it)      = sum(spks(it,:) > (Task.StimOnset(it) + preVisWind(1)) & spks(it,:) < Task.StimOnset(it) + preVisWind(2))/(diff(preVisWind)/1000);
    visTrans(it)    = sum(spks(it,:) > (Task.StimOnset(it) + visTransWind(1)) & spks(it,:) < Task.StimOnset(it) + visTransWind(2))/(diff(visTransWind)/1000);
    visSust(it)     = sum(spks(it,:) > (Task.StimOnset(it) + visSustWind(1)) & spks(it,:) < Task.StimOnset(it) + visSustWind(2))/(diff(visSustWind)/1000);
    prePreSacc(it)  = sum(spks(it,:) > (Task.SRT(it) + Task.GoCue(it) + prePreSaccWind(1)) & spks(it,:) < Task.SRT(it) + Task.GoCue(it) + prePreSaccWind(2))/(diff(prePreSaccWind)/1000);
    preSacc(it)     = sum(spks(it,:) > (Task.SRT(it) + Task.GoCue(it) + preSaccWind(1)) & spks(it,:) < (Task.SRT(it) + Task.GoCue(it) + preSaccWind(2)))/(diff(preSaccWind)/1000);
    postSacc(it)    = sum(spks(it,:) > (Task.SRT(it) + Task.GoCue(it) + postSaccWind(1)) & spks(it,:) < Task.SRT(it) + Task.GoCue(it) + postSaccWind(2))/(diff(postSaccWind)/1000);
    postTone(it)    = sum(spks(it,:) > (Task.RewardTone(it) + postToneWind(1)) & spks(it,:) > Task.RewardTone(it) + postToneWind(2))/(diff(postToneWind)/1000);
    postRewd(it)    = sum(spks(it,:) > (Task.Reward(it) + postRewdWind(1)) & spks(it,:) < (Task.Reward(it) + postRewdWind(2)))/(diff(postRewdWind)/1000);
end

[errBL, preErr, postErr] = deal(nan(1,size(spks,1)));
for ii = 1:length(errTr),
    it = errTr(ii);
    errBL(it)     = sum(spks(it,:) > (Task.StimOnset(it) + preVisWind(1)) & spks(it,:) < (Task.StimOnset(it) + preVisWind(2)))/(diff(preVisWind)/1000);
    preErr(it)    = sum(spks(it,:) > (Task.StimOnset(it) + preErrWind(1)) & spks(it,:) < (Task.StimOnset(it) + preErrWind(2)))/(diff(preErrWind)/1000);
    postErr(it)   = sum(spks(it,:) > (Task.StimOnset(it) + postErrWind(1)) & spks(it,:) < (Task.StimOnset(it) + postErrWind(2)))/(diff(postErrWind)/1000);
end

% Get visuomotor index (   (visRate - baseRate)/(movRate - baseRate)    )
%      Cohen et al 2009 Definition
visRate = nanmean(visTrans); movRate = nanmean(preSacc); baseRate = nanmean(preVis);

vmi = (visRate - baseRate)/(movRate - baseRate);

%% Do comparisons
% Baseline = PreVis
% bl vs visTrans
% bl vs visSust
% bl vs preSacc
% preSacc vs prePreSacc
% bl vs postSacc
% visTrans vs postSacc
% bl vs postTone
% bl vs postRewd
% 
% errBL vs preErr
% errBL vs postErr

pVals(1) = signrank(preVis,visTrans);
pVals(2) = signrank(preVis,visSust);
pVals(3) = signrank(preVis,preSacc);
pVals(4) = signrank(prePreSacc,preSacc);
pVals(5) = signrank(preVis,postSacc);
pVals(6) = signrank(visTrans,postSacc);
pVals(7) = signrank(preVis,postTone);
pVals(8) = signrank(preVis,postRewd);
pVals(9) = signrank(errBL,preErr);
pVals(10) = signrank(errBL,postErr);

compTypes = {'vTrans','vSust','pSacc','VSvsPreSacc','postSacc','VTvsPostSacc','rwdTone','Rwd','preErr','Err'};

%% I used the following line for correction for multiple comparisons but commented it out because I didn't write the pAdjust function I use
pOut = pAdjust(pVals);

% So instead this will export the uncorrected p values - might be
% problematic for vis/mov identifications?
%pOut = pVals;

% Let's define visual index as mean rate 50-150ms post stimulus/baseline
% and motor index = mean rate 100ms pre-saccade/baseline
zVis = (nanmean(visTrans)-nanmean(preVis))./(nanstd(preVis)./sqrt(length(preVis)));
zMov = (nanmean(preSacc)-nanmean(prePreSacc))./(nanstd(prePreSacc)./sqrt(length(prePreSacc)));

vmi = [zVis,zMov];

spkRate = [preVis;visTrans;visSust;preSacc;postSacc;postTone;postRewd; errBL; preErr; postErr];
rateCols = {'PreVis','VisTrans','VisSust','PreSacc','PostSacc','PostTone','PostRewd','ErrBaseline','PreErr','PostErr'};

isVis = 0; isMov = 0;
if pOut(1) < .05 || pOut(2) < .05, isVis = 1; end
if pOut(3) < .05 && pOut(4) < .05, isMov = 1; end
% if pOut(3) < .05 && pOut(4) > .05
%     keyboard
% end
if isVis && isMov,
    type = 'vismov';
elseif isVis && ~isMov,
    type = 'vis';
elseif ~isVis && isMov,
    type = 'mov';
else
    type = 'none';
end


end
