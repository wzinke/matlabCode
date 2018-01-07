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

% load([tebaMount,'Users/Kaleb/proAntiProcessed/Darwin-170830-141157/Behav.mat');

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
% nBeta = 2;
nPis = 5;
% nBeta = (length(params)-nPis-1)/2;
nBeta = (length(params)-nPis)/2;
antiInd = 1;
proInd = 2;
% betas = params(((end-(nBeta-1)):end)-1);
% pis = params(1:(end-(nBeta+1)));
pis = params(1:nPis);
betas = params(nPis+(1:(nBeta*2)));
% lamb = params(end);
lamb = 1;

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
    
%     respChars = stimChars(:,(end-(nBeta-1)):end);
%     stimChars = stimChars(:,1:(end-(nBeta)));
    respChars = stimChars(:,(nPis+1):end);
    stimChars = stimChars(:,1:nPis);
    

    % Use stimulus characteristics and pis to get a probability of
    % selecting a stimulus
    pSel = stimChars*pis;
    pSel = pSel./nansum(pSel);
    
    % Probability of looking at a stimulus (pLook(x)) =
    % sum(pLook(x|y)*pSel(y)) for all ys
    % In TVA terms: sum(eta(i,x)*beta(x,y)*eta(i,y)
    % Rather, evidence that i should be associated with a given y's beta(x) *
    % that beta * evidence that i is y
    
    % Let's set up myFactors
    stimInds = 1:length(pSel);
    myFactors = nan(length(pSel),nBeta,length(pSel));
    for iSel = 1:length(pSel)
        myFactors(:,1:2,iSel) = respChars(:,1:2);
        for i = 0:(length(pSel)/2)
            myFactors(:,i+3,iSel) = ((circshift(stimInds,i)-iSel)==0)|((circshift(stimInds,-i)-iSel)==0);
            myFactors(myFactors(:,i+3,iSel)==0,i+3,iSel) = .1;
        end
        myFactors(:,[6,7],iSel) = respChars;
    end
    
    % Structure betas: 1 is horizontal stimulus, 2 is vertical stimulus
    myBetas(1,:) = betas(1:nBeta);
    myBetas(2,:) = betas((1:nBeta)+nBeta);
    
    % Selected stimulus loop
    pLookTmp = nan(length(pSel),length(pSel));
    for iSel = 1:length(pSel)
        % All stimulus loop
        stimLook = zeros(length(pSel),1);
        for iStim = 1:length(pSel)
            % Beta loop: vertical selected or horizontal selected
            for iBeta = 1:2
                % Factor loop
                for iFactor = 1:size(myBetas,2)
                    % iStim is iFactor*myBetas(iBeta,iFactor)*iSel is V or
                    % H
                    stimLook(iStim) = stimLook(iStim)+myFactors(iStim,iFactor,iSel)*myBetas(iBeta,iFactor)*respChars(iSel,iBeta);
                end
            end
        end
        pLookTmp(:,iSel) = stimLook./nansum(stimLook);
    end
%     pLookTmp = pLookTmp./nansum(pLookTmp(:));                
    pLook = pLookTmp*pSel;
    pLook = pLook./nansum(pLook);
    
    err(it) = sum(abs(pLook-stimChosen));
end
prob = nansum(err);

