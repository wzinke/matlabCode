function [prob, pLook] = klTVA(params,Task)


% What we'll do is have an mxn matrix of stimulus characteristics where m
% is the set size and n is the number of characteristics. Thus, if we
% multiply by column vector pis, we will get an mx1 vector of stimulus
% weights wx
reverse = 0;
if sum(size(params) > 1) > 1
    error('Params needs to be a column vector, not a matrix');
elseif size(params,1) == 1
    params = params';
    reverse = 1;
end

% load('/mnt/teba/Users/Kaleb/proAntiProcessed/Darwin-170830-141157/Behav.mat');

stimHs = [.5 1 2];
stimVs = [2 1 .5];

%% Set up stimChars
isPA = strcmpi(Task.TaskType,'Pro-Anti');
isGood = Task.Abort == 0;
Task = klCutTask(Task,isPA & isGood);
nTrs = length(Task.trStarts);
err = nan(nTrs,1);
% pisRaw = pis;

%% Set up pis and betas from "params"
nBetaParams = 1;
nBeta = 2;
antiInd = 1;
proInd = 2;
betaParams = params((end-(nBetaParams-1)):end);
piParams = params(1:(end-(nBetaParams)));
betas = [betaParams(1),1-betaParams(1)]';
pis = [piParams(1),1-piParams(1)]';


for it = 1:nTrs
    % Stimchars = mxn m= set size, n = red, green, vertical, horizontal,
    % square
    if isnan(Task.EndStimInd(it))
        continue
    else
        stimChosen = zeros(Task.SetSize(it),1);
        stimChosen(Task.EndStimInd(it)) = 1;
    end    
    
    stimChars = nan(Task.SetSize(it),4);
    for is = 1:Task.SetSize(it)
        if Task.SingletonColor(it) == 0 && Task.StimLoc(it,is)==Task.TargetLoc(it)
            stimChars(is,1) = 1;
            stimChars(is,2) = 0;
        else
            stimChars(is,1) = 0;
            stimChars(is,2) = 1;
        end
        stimChars(is,3) = normcdf(log(stimHs(Task.StimDiff(it,is)+1)/stimVs(Task.StimDiff(it,is)+1)));
        stimChars(is,4) = normcdf(log(stimVs(Task.StimDiff(it,is)+1)/stimHs(Task.StimDiff(it,is)+1)));
%         stimChars(is,5) = 1-abs(normcdf(log(stimHs(Task.StimDiff(it,is)+1)/stimVs(Task.StimDiff(it,is)+1)))-.5)*2;
    end
    
    
%     stimChars = [0 1 .5 .5 0;...
%        1 0 .8 .2 0;...
%        0 1 .5 .5 0;...
%        0 1 .5 .5 0];
    
%     stimChars = [0 1 0 1 0;...
%         1 0 1 0 0;...
%         0 1 0 1 0;...
%         0 1 0 1 0];
% 
%     stimChars = [0 1 normcdf(log(.5/2),0,1) normcdf(log(2/.5),0,1) (1-abs(normcdf(log(2/.5))-.5)*2);...
%         1 0 normcdf(log(2/.5),0,1) normcdf(log(.5/2),0,1) (1-abs(normcdf(log(2/.5))-.5)*2);...
%         0 1 normcdf(log(.5/2),0,1) normcdf(log(2/.5),0,1) (1-abs(normcdf(log(2/.5))-.5)*2);...
%         0 1 normcdf(log(.5/2),0,1) normcdf(log(2/.5),0,1) (1-abs(normcdf(log(2/.5))-.5)*2)];
%     
    


    respChars = stimChars(:,(end-(nBeta-1)):end);
    stimChars = stimChars(:,1:(end-(nBeta)));

    %% Note this will be for one trial, it seems
    wx = stimChars*pis;

    % Get vx = eta * beta * weights
    vx = respChars.*repmat(betas,1,length(wx))'.*repmat(wx./sum(wx),1,nBeta);

    % Make vx sum to 1
    px = vx./sum(vx(:));

    % Get pLook(stim) = px(stim,pro) + px(opp,anti)
    antiLocs = mod((0:(length(wx)-1))+(length(wx)/2),length(wx))+1;
    pLook = px(:,proInd)+px(antiLocs,antiInd);
    err(it) = sum(abs(pLook-stimChosen));
end
prob = nansum(err);
