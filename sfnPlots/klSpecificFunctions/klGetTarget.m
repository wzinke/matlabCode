function   [pAdj, pNames] = klGetTarget(spks,Task,varargin)

% Decode varargin
if nargin > 4
    % Check that there are enough input arguments, if optional flags are
    % set
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

%% Define intervals
preTrialWind    = [-100, 0];
preVisWind      = [-100, 0];
visTransWind    = [0, 100];
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
[tLoc dLoc sing] = deal(nan(1,size(spks,1)));

corrTr = find(Task.Correct);
errTr  = find(~isnan(Task.ErrorTone));

%% Get Rates
for ii = 1:length(corrTr),
    it = corrTr(ii); 
    preTrial(it)    = sum(spks(it,:) > preTrialWind(1) & spks(it,:) < preTrialWind(2))/(diff(preTrialWind)/1000);
    preVis(it)      = sum(spks(it,:) > (Task.StimOnset(it) + preVisWind(1)) & spks(it,:) < Task.StimOnset(it) + preVisWind(2))/(diff(preVisWind)/1000);
    visTrans(it)    = sum(spks(it,:) > (Task.StimOnset(it) + visTransWind(1)) & spks(it,:) < Task.StimOnset(it) + visTransWind(2))/(diff(visTransWind)/1000);
    visSust(it)     = sum(spks(it,:) > (Task.StimOnset(it) + visSustWind(1)) & spks(it,:) < Task.StimOnset(it) + visSustWind(2))/(diff(visSustWind)/1000);
    prePreSacc(it)  = sum(spks(it,:) > (Task.SRT(it) + prePreSaccWind(1)) & spks(it,:) < Task.SRT(it) + prePreSaccWind(2))/(diff(prePreSaccWind)/1000);
    preSacc(it)     = sum(spks(it,:) > (Task.SRT(it) + preSaccWind(1)) & spks(it,:) < (Task.SRT(it) + preSaccWind(2)))/(diff(preSaccWind)/1000);
    postSacc(it)    = sum(spks(it,:) > (Task.SRT(it) + postSaccWind(1)) & spks(it,:) < Task.SRT(it) + postSaccWind(2))/(diff(postSaccWind)/1000);
    postTone(it)    = sum(spks(it,:) > (Task.RewardTone(it) + postToneWind(1)) & spks(it,:) > Task.RewardTone(it) + postToneWind(2))/(diff(postToneWind)/1000);
    postRewd(it)    = sum(spks(it,:) > (Task.Reward(it) + postRewdWind(1)) & spks(it,:) < (Task.Reward(it) + postRewdWind(2)))/(diff(postRewdWind)/1000);
    tLoc(it) = Task.TargetLoc(it);
    if isfield(Task,'DistLoc'), dLoc(it) = Task.DistLoc(it); end
    if isfield(Task,'Singleton'), sing(it) = Task.Singleton(it); end
end
pTarg(1) = anovan(visTrans,{tLoc},'varnames',{'TargetLoc'},'display','off');
pTarg(2) = anovan(visSust,{tLoc},'varnames',{'TargetLoc'},'display','off');
pTarg(3) = anovan(preSacc,{tLoc},'varnames',{'TargetLoc'},'display','off');
pTarg(4) = anovan(postSacc,{tLoc},'varnames',{'TargetLoc'},'display','off');

allPs(1,1:4) = pTarg(1:4);
pNames = {'VT-Targ','VS-Targ','PreS-Targ','PostS-Targ'};

uniqueLocs = unique(tLoc(isfinite(tLoc)));
for il = 1:length(uniqueLocs)
    trialsType{il,1} = tLoc == uniqueLocs(il);
    trialsType{il,2} = dLoc == uniqueLocs(il) & sing == 1;
    trialsType{il,3} = dLoc == uniqueLocs(il) & sing == 0;
    
    aovType = [];
    aovTrans = [];
    aovSust = [];
    aovPreSacc = [];
    aovPostSacc = [];
    
    for it = 1:3
        aovTrans = cat(2,aovTrans,visTrans(trialsType{il,it}));
        aovSust = cat(2,aovSust,visSust(trialsType{il,it}));
        aovPreSacc = cat(2,aovPreSacc,preSacc(trialsType{il,it}));
        aovPostSacc = cat(2,aovPostSacc,postSacc(trialsType{il,it}));
        aovType = cat(2,aovType,repmat(it,1,sum(trialsType{il,it})));
    end
    pTransType(il,1) = anovan(aovTrans,{aovType},'display','off');
    pTransType(il,2) = anovan(aovSust,{aovType},'display','off');
    pTransType(il,3) = anovan(aovPreSacc,{aovType},'display','off');
    pTransType(il,4) = anovan(aovPostSacc,{aovType},'display','off');
    
    allPs = cat(2,allPs,pTransType(il,:));
    pNames = cat(2,pNames,{sprintf('VT-Loc %d', uniqueLocs(il)),sprintf('VS-Loc %d', uniqueLocs(il)),sprintf('PreS-Loc %d', uniqueLocs(il)),sprintf('PostS-Loc %d', uniqueLocs(il))});
    
end
pAdj = pAdjust(allPs);

%keyboard