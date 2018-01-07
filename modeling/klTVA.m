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

% load([tebaMount,'Users/Kaleb/proAntiProcessed/Darwin-170830-141157/Behav.mat']);

stimHs = [.5 1 2];
stimVs = [2 1 .5];

%% Set up stimChars
if isfield(Task,'TaskType')
    isPA = strcmpi(Task.TaskType,'Pro-Anti');
else
    isPA = ones(length(Task.Abort),1);
end
isGood = Task.Abort == 0;
Task = klCutTask(Task,isPA & isGood);
nTrs = length(Task.Correct);
err = nan(nTrs,1);
% pisRaw = pis;

%% Set up pis and betas from "params"
nBeta = 2;
antiInd = 1;
proInd = 2;
betas = params(((end-(nBeta-1)):end)-1);
pis = params(1:(end-(nBeta+1)));
lamb = params(end);

for it = 1:nTrs
    % Stimchars = mxn m= set size, n = red, green, vertical, horizontal,
    % square
    if isnan(Task.EndStimInd(it))
        continue
    else
        stimChosen = zeros(Task.SetSize(it),1);
        stimChosen(Task.EndStimInd(it)) = 1;
    end    
    
    oppLocs = mod((0:(Task.SetSize(it)-1))+(Task.SetSize(it)/2),Task.SetSize(it))+1;

    stimChars = nan(Task.SetSize(it),4);
    for is = 1:Task.SetSize(it)
        % 1 and 2: Color
        if Task.SingletonColor(it) == 0 && Task.StimLoc(it,is)==Task.TargetLoc(it)
            stimChars(is,1) = 1;
            stimChars(is,2) = .1;
        else
            stimChars(is,1) = .1;
            stimChars(is,2) = 1;
        end
        % 3-5: Aspect ratios        
        stimChars(is,3) = normcdf(log(stimHs(Task.StimDiff(it,is)+1)/stimVs(Task.StimDiff(it,is)+1)),0,lamb);
        stimChars(is,4) = normcdf(log(stimVs(Task.StimDiff(it,is)+1)/stimHs(Task.StimDiff(it,is)+1)),0,lamb);
        stimChars(is,5) = 1-abs(normcdf(log(stimHs(Task.StimDiff(it,is)+1)/stimVs(Task.StimDiff(it,is)+1)),0,lamb)-.5)*2;
        % 6-7: Horizontal (anti) and Vertical (pro)
        stimChars(is,6) = normcdf(log(stimHs(Task.StimDiff(it,is)+1)/stimVs(Task.StimDiff(it,is)+1)),0,lamb);
        stimChars(is,7) = normcdf(log(stimVs(Task.StimDiff(it,is)+1)/stimHs(Task.StimDiff(it,is)+1)),0,lamb);
        %8:7+SetSize: Location
%         stimChars(is,7+(1:Task.SetSize(is))) = 0; stimChars(is,7+is) = 1;
    end
    
    respChars = stimChars(:,(end-(nBeta-1)):end);
    stimChars = stimChars(:,1:(end-(nBeta)));

    %% Note this will be for one trial, it seems
    wx = stimChars*pis;

    vx(:,1) = respChars(:,1).*betas(1).*(wx(oppLocs)./sum(wx));
    vx(:,2) = respChars(:,2).*betas(2).*(wx./sum(wx));
%     vx = vx.*repmat(wx./sum(wx),1,nBeta);
%     vx = respChars.*repmat(betas,1,length(wx))'.*repmat(wx./sum(wx),1,nBeta);

    % Make vx sum to 1
    px = vx./sum(vx(:));

    % Get pLook(stim) = px(stim,pro) + px(opp,anti)
%     antiLocs = mod((0:(length(wx)-1))+(length(wx)/2),length(wx))+1;
%     pLook = px(:,proInd)+px(antiLocs,antiInd);
    pLook = nansum(px,2);
    
% % %     % Get vx = eta * beta * weights
% % %     vx = respChars.*repmat(betas,1,length(wx))'.*repmat(wx./sum(wx),1,nBeta);
% % % 
% % %     % Make vx sum to 1
% % %     px = vx./sum(vx(:));
% % % 
% % %     % Get pLook(stim) = px(stim,pro) + px(opp,anti)
% % %     antiLocs = mod((0:(length(wx)-1))+(length(wx)/2),length(wx))+1;
% % %     pLook = px(:,proInd)+px(antiLocs,antiInd);
    err(it) = sum(abs(pLook-stimChosen));
end
prob = nansum(err);

