function [grpSDF,grpSDFTimes,normMat,normMatID] = klGetSDFs(spiketimes,Task)

% Set defaults
blWind = -200:-100;
vWind = -100:400;
mWind = -400:100;
zDim = [1,2];
zType = 'baseline';
whichTask = 'MG';
doRF = 0;
groupVar = 'none';
clip = 1;

% Decode varargin

% Set saccade times
saccTime = Task.GoCue+Task.SRT;

% Make the vis-aligned SDF
checkGood = spiketimes(ismember(Task.TaskType,'MG'),:);
if sum(isnan(checkGood(:))) == numel(checkGood),
    normMatID = {'blMean','blStd','trMean','trStd','vAbsMax','mAbsMax'};
    normMat = nan(1,6);
    [grpSDF{1,1:3}] = deal(nan);
    [grpSDFTimes{1,1:3}] = deal(nan);
    return;
end

[tmpSDF,tmpTimes] = klSpkRatev2(spiketimes(ismember(Task.TaskType,whichTask),:),'-q',1);

% Get baseline rate
bl = nanmean(nanmean(tmpSDF(:,ismember(tmpTimes,blWind)),zDim(1)),zDim(2));

% Get the best-guess for RF in case we want a "maximally driven" SDF
[rf ant1f]= getRFv2(spiketimes,Task,'-b',bl);
vRF = rf{1}; mRF = rf{2};
vARF = ant1f(1); mARF = ant1f(2);
if doRF,
    vCrit = ismember(Task.TargetLoc,vRF);
    mCrit = ismember(Task.TargetLoc,mRF);
    vaCrit = ismember(Task.TargetLoc,vARF);
    maCrit = ismember(Task.TargetLoc,mARF);
else
    vCrit = ismember(Task.TaskType,whichTask) & Task.Correct==1;
    mCrit = ismember(Task.TaskType,whichTask) & Task.Correct==1;
end

% Get a visually aligned SDF from the appropriate trials
[rawVisSDF,grpSDFTimes{1,3}] = klSpkRatev2(spiketimes(vCrit,:),'-q',1);
% Get saccade-aligned SDF
[mSDF, grpSDFTimes{1,2}] = klSpkRatev2(spiketimes(mCrit,:)-repmat(Task.SRT(mCrit,:)+Task.GoCue(mCrit,:),1,size(spiketimes,2)),'-q',1);

% Cut out post-saccade spikes if we want
if clip,
    spiketimes(spiketimes > repmat(saccTime,1,size(spiketimes,2))) = nan;
end

% Get visually aligned SDF with the post-saccade spikes removed
[vSDF, grpSDFTimes{1,1}] = klSpkRatev2(spiketimes(vCrit,:),'-q',1);

% Here is the baseline mean, SD, visual max, and movement max
normMatID = {'blMean','blStd','trMean','trStd','vAbsMax','mAbsMax'};
normMat(1,1) = nanmean(nanmean(vSDF(:,ismember(grpSDFTimes{1,1},blWind)),zDim(1)),zDim(2));
normMat(1,2) = nanstd(nanmean(rawVisSDF(:,ismember(grpSDFTimes{1,1},blWind)),zDim(1)),[],zDim(2));
normMat(1,3) = nanmean(nanmean(rawVisSDF,zDim(1)),zDim(2));
normMat(1,4) = nanstd(nanmean(rawVisSDF,zDim(1)),[],zDim(2));
normMat(1,5) = nanmax(abs(nanmean(vSDF(:,ismember(grpSDFTimes{1,1},vWind)),1)-normMat(1,1)));
normMat(1,6) = nanmax(abs(nanmean(mSDF(:,ismember(grpSDFTimes{1,2},mWind)),1)-normMat(1,1)));

% If we want to organize by some trial grouping, do that here
switch groupVar
    case 'none'
        vTrGroup = ones(size(vSDF,1),1);
        mTrGroup = ones(size(mSDF,1),1);
    case {'loc'}
        vTrGroup = Task.TargetLoc;
        mTrGroup = Task.TargetLoc;
    case {'targ','target','tst'}
        vTrGroup = ismember(Task.TargetLoc,vRF);
        mTrGroup = ismember(Task.TargetLoc,mRF);
end

% Place the mean SDF as organized by the trial grouping
uGroup = unique(vTrGroup);
for ig = 1:length(uGroup)
    grpSDF{1,1}(ig,:) = nanmean(vSDF(vTrGroup == uGroup(ig),:),1);
    grpSDF{1,3}(ig,:) = nanmean(rawVisSDF(vTrGroup == uGroup(ig),:),1);
    grpStd{1,1}(ig,:) = nanstd(vSDF(vTrGroup == uGroup(ig),:),[],1)./sqrt(sum(isfinite(vSDF(vTrGroup == uGroup(ig),:)),1));
    grpStd{1,3}(ig,:) = nanstd(rawVisSDF(vTrGroup == uGroup(ig),:),[],1)./sqrt(sum(isfinite(rawVisSDF(vTrGroup == uGroup(ig),:)),1));
end
uGroup = unique(mTrGroup);
for ig = 1:length(uGroup)
    grpSDF{1,2}(ig,:) = nanmean(mSDF(mTrGroup == uGroup(ig),:),1);
    grpStd{1,2}(ig,:) = nanstd(mSDF(mTrGroup == uGroup(ig),:),[],1)./sqrt(sum(isfinite(mSDF(mTrGroup == uGroup(ig),:)),1));
end
    