function klViewSearch(allSpikes,Task)

% Set constants
tLocs = 0:90:270;
tSubInds = [6,2,4,8];
vWind = [-300, 400];
mWind = [-400,300];
clip = 0;

% First, get MG trials
if ~any(ismember(Task.TaskType,{'Cap','Search'})),
    warning('No Search trials found...');
    return
else
    capTrials = ismember(Task.TaskType,{'Cap','Search'});
    spikes = allSpikes(capTrials,:);
    
    tFields = fieldnames(Task);
    goodLen = size(Task.Correct,1);
    for iField = 1:length(tFields),
        if size(Task.(tFields{iField}),1) == goodLen,
            capTask.(tFields{iField}) = Task.(tFields{iField})(capTrials,:);
        else
            capTask.(tFields{iField}) = Task.(tFields{iField});
        end
    end
end

% Make SDFs
[mSDF,mTimes] = klSpkRatev3(spikes-repmat(capTask.SRT+capTask.GoCue,1,size(spikes,2)),'-q',1);

% Cut post-saccade spikes
if clip,
    cutSpks = spikes;
    cutSpks(cutSpks >= repmat(capTask.SRT+capTask.GoCue,1,size(cutSpks,2))) = nan;
else
    cutSpks = spikes;
end
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
    pltMeanStd(vTimes,nanmean(vSDF(capTask.TargetLoc==tLocs(il),:),1),nanstd(vSDF(capTask.TargetLoc==tLocs(il),:),[],1)./sqrt(sum(capTask.TargetLoc==tLocs(il))),'color','r');
    pltMeanStd(vTimes,nanmean(vSDF(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==1,:),1),nanstd(vSDF(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==1,:),[],1)./sqrt(sum(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==1)),'color','g');
    pltMeanStd(vTimes,nanmean(vSDF(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==0,:),1),nanstd(vSDF(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==0,:),[],1)./sqrt(sum(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==0)),'color','b');
    vline(0);
    set(gca,'XLim',vWind);
    
    % Movement response:
    figure(mFig); subplot(3,3,tSubInds(il));
    pltMeanStd(mTimes,nanmean(mSDF(capTask.TargetLoc==tLocs(il),:),1),nanstd(mSDF(capTask.TargetLoc==tLocs(il),:),[],1)./sqrt(sum(capTask.TargetLoc==tLocs(il))),'color','r');
    pltMeanStd(mTimes,nanmean(mSDF(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==1,:),1),nanstd(mSDF(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==1,:),[],1)./sqrt(sum(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==1)),'color','g');
    pltMeanStd(mTimes,nanmean(mSDF(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==0,:),1),nanstd(mSDF(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==0,:),[],1)./sqrt(sum(capTask.TargetLoc~=tLocs(il) & capTask.Singleton==0)),'color','b');
    vline(0);
    set(gca,'XLim',mWind);
end