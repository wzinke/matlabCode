function klViewMG(allSpikes,Task)

% Set constants
tLocs = 0:45:315;
tSubInds = [6,3,2,1,4,7,8,9];
vWind = [-300, 700];
mWind = [-700,300];

% First, get MG trials
if ~any(ismember(Task.TaskType,'MG')),
    warning('No MG trials found...');
    return
else
    mgTrials = strcmp(Task.TaskType,'MG');
    spikes = allSpikes(mgTrials,:);
    
    tFields = fieldnames(Task);
    goodLen = size(Task.Correct,1);
    for iField = 1:length(tFields),
        if size(Task.(tFields{iField}),1) == goodLen,
            mgTask.(tFields{iField}) = Task.(tFields{iField})(mgTrials,:);
        else
            mgTask.(tFields{iField}) = Task.(tFields{iField});
        end
    end
end

% Make SDFs
[mSDF,mTimes] = klSpkRatev3(spikes-repmat(mgTask.SRT+mgTask.GoCue,1,size(spikes,2)),'-q',1);

% Cut post-saccade spikes
cutSpks = spikes;
cutSpks(cutSpks >= repmat(mgTask.SRT+mgTask.GoCue,1,size(cutSpks,2))) = nan;
[vSDF,vTimes] = klSpkRatev3(cutSpks,'-q',1);

% Open figures:
if ~exist('vFig','var'),
    vFig = figure();
end
if ~exist('mFig','var'),
    mFig = figure();
end

% Loop through target locations
for il = 1:length(tLocs),
    % Visual response:
    figure(vFig); subplot(3,3,tSubInds(il));
    pltMeanStd(vTimes,nanmean(vSDF(mgTask.TargetLoc==tLocs(il),:),1),nanstd(vSDF(mgTask.TargetLoc==tLocs(il),:),[],1)./sqrt(sum(mgTask.TargetLoc==tLocs(il))));
    vline(0);
    set(gca,'XLim',vWind);
    
    % Movement response:
    figure(mFig); subplot(3,3,tSubInds(il));
    pltMeanStd(mTimes,nanmean(mSDF(mgTask.TargetLoc==tLocs(il),:),1),nanstd(mSDF(mgTask.TargetLoc==tLocs(il),:),[],1)./sqrt(sum(mgTask.TargetLoc==tLocs(il))));
    vline(0);
    set(gca,'XLim',mWind);
end