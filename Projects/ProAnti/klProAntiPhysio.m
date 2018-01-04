function klProAntiPhysio(sessName,varargin)

baseDir = '/mnt/teba/Users/Kaleb/proAntiProcessed/';
doSave = 1;
watch = 0;
fresh = 1;
channel = 0;
pause = 0;

varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-s'}
            doSave = varargin{varStrInd(iv)+1};
        case {'-f'}
            fresh = varargin{varStrInd(iv)+1};
        case {'-c'}
            channel = varargin{varStrInd(iv)+1};
        case {'-p'}
            pause = varargin{varStrInd(iv)+1};
    end
end

% Load in behavior and eye positions
% sessDir = chanStruct.folder;
sessDir = [baseDir,filesep,sessName];
% [~,sessName] = fileparts(sessDir);
load([sessDir, '/Behav.mat']);
% load([sessDir, '/Eyes.mat']);

% Find unit .mat files
% chanMats = dir([sessDir,'/chan*.mat']);
allUnits = dir([sessDir,filesep,'Chan*']);
unitNames = {allUnits.name};
if ~isnumeric(channel) || channel ~=0
    analyzeChannels = unitNames(ismember(unitNames,channel));
else
    analyzeChannels = unitNames;
end

% Set window for visual RF determination
vWind = 50:150;
vXwind = [-200,400];
mXwind = [-400,200];


% Other helpful stuff
congStyle = {'-','--',':'};

% Make some logicals
isPA = strcmpi(Task.TaskType,'Pro-Anti');
isPro = strcmp(Task.TrialType,'pro');
isAnti = strcmp(Task.TrialType,'anti');
isCong = Task.Congruent == 0;
isIncong = Task.Congruent == 1;
isCatch = Task.Congruent == 2;
congMat = [isCong,isIncong,isCatch];
isCorrect = Task.Correct == 1;
isGood = Task.Abort == 0;

shapeSingMat = klUniqueMat(Task.StimDiff);

endDiff = nan(length(Task.EndStimInd),1);
for it = 1:length(Task.EndStimInd)
    if ~isnan(Task.EndStimInd)
        endDiff(it) = Task.StimDiff(it,Task.EndStimInd(it));
    end
end

% Get unique target locations
proLocs = nunique(Task.TargetLoc(isPA));
antiLocs = mod(proLocs+180,360);
uEccs = nunique(Task.Eccentricity(isPA));
if length(proLocs)==4
    spInds = [6,2,4,8];
end

keepVars = whos; keepVars = {keepVars.name}; keepVars{end+1} = 'keepVars'; keepVars{end+1} = 'clrVars'; keepVars{end+1} = 'ic';

% Start unit loop
for ic = 1:length(analyzeChannels)
    close all;
    currVars = whos; currVars = {currVars.name}; 
    clrVars = currVars(~ismember(currVars,keepVars)); 
    for iv = 1:length(clrVars), clear(clrVars{iv}); end
    
    fprintf('Working on unit %d (of %d)...\n',ic,length(analyzeChannels));
    saveDir = [sessDir,filesep,analyzeChannels{ic}];
    if doSave && length(dir(sprintf('%s/PA*.png',saveDir))) == 7 && ~fresh
        continue
    end

    load([sessDir,filesep,analyzeChannels{ic},filesep,lower(analyzeChannels{ic}(1)),analyzeChannels{ic}(2:end),'.mat']);
    spkMat = klPlaceEvents(Task,spkTimes);

    if size(spkMat,2) < 5
        continue
    end

    [vSDF,vTimes] = klSpkRatev2(spkMat,'-q',1,'-k','psp');
    [mSDF,mTimes] = klSpkRatev2(spkMat-repmat(Task.SRT,1,size(spkMat,2)),'-q',1,'-k','psp');

    % Figure 1 will plot the visual SDF in all locations for pro and anti
    % trials (correct trials where the saccade is made TOWARD that
    % location)
    figure(1);
    set(1,'renderer','zbuffer');
    for il = 1:length(proLocs)
        % Good pro trials
        proTrials(:,il)   = isPA & isPro & isCorrect & Task.TargetLoc == proLocs(il);
        antiTrials(:,il)  = isPA & isAnti & isCorrect & Task.TargetLoc == proLocs(il);

        % Output
        mnProSDFV(il,:) = nanmean(vSDF(proTrials(:,il),:),1);
        mnAntiSDFV(il,:) = nanmean(vSDF(antiTrials(:,il),:),1);
        mnProSDFM(il,:) = nanmean(mSDF(proTrials(:,il),:),1);
        mnAntiSDFM(il,:) = nanmean(mSDF(antiTrials(:,il),:),1);

        sdProSDFV(il,:) = nanstd(vSDF(proTrials(:,il),:),[],1);
        sdAntiSDFV(il,:) = nanstd(vSDF(antiTrials(:,il),:),[],1);

        sp(il) = subplot(3,3,spInds(il)); hold on;

        plot(vTimes,mnProSDFV(il,:),'linewidth',3,'color','k');
        plot(vTimes,mnAntiSDFV(il,:),'linewidth',1,'color','k');
    end
    set(sp,'XLim',vXwind);
    linkaxes(sp,'y');
    yLims = get(sp(1),'YLim');
    for il = 1:length(proLocs)
        subplot(sp(il));
        plot(sort(Task.SRT(proTrials(:,il) & ~isnan(Task.SRT))),(1:sum(proTrials(:,il) & ~isnan(Task.SRT))).*yLims(2)./sum(proTrials(:,il) & ~isnan(Task.SRT)),'color','g');
        plot(sort(Task.SRT(antiTrials(:,il) & ~isnan(Task.SRT))),(1:sum(antiTrials(:,il) & ~isnan(Task.SRT))).*yLims(2)./sum(antiTrials(:,il) & ~isnan(Task.SRT)),'color','r');
        set(gca,'children',flipud(get(gca,'children')));
    end
    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(1,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(1,sprintf('%s/PA-VisResponses.png',saveDir));
        close(1);
    end

    % Figure 2 is the same as above, but aligned on saccade
    figure(2);
    set(2,'renderer','zbuffer');
    for il = 1:length(proLocs)
        sp(il) = subplot(3,3,spInds(il)); hold on;
        plot(mTimes,mnProSDFM(il,:),'linewidth',3,'color','k');
        plot(mTimes,mnAntiSDFM(il,:),'linewidth',1,'color','k');
    end
    set(sp,'XLim',mXwind);
    linkaxes(sp,'y');
    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(2,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(2,sprintf('%s/PA-MovResponses.png',saveDir));
        close(2);
    end

    % Calculate RF from pro trials
    [~,maxInd] = max(nanmean(mnProSDFV(:,ismember(vTimes,vWind)),2));
    oppLoc = mod(maxInd+(length(proLocs)/2),length(proLocs));
    if oppLoc == 0, oppLoc = length(proLocs); end
    vRF = proLocs(maxInd);
    aRF = proLocs(oppLoc);

    % Figure 3 will plot the activity when the eye movement is made toward
    % the RF across the different congruency conditions
    figure(3); hold on;
    set(3,'renderer','zbuffer');
    for ii = 1:3
        plot(vTimes,nanmean(vSDF(proTrials(:,maxInd) & congMat(:,ii),:),1),'color','g','linewidth',1,'linestyle',congStyle{ii});
        plot(vTimes,nanmean(vSDF(antiTrials(:,oppLoc) & congMat(:,ii),:),1),'color','r','linewidth',1,'linestyle',congStyle{ii});
    end
    set(gca,'XLim',vXwind);
    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(3,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(3,sprintf('%s/PA-Cong-SaccDir.png',saveDir));
        close(3);
    end

    % Figure 4 will plot the activity when the stimulus is in
    % the RF across the different congruency conditions
    figure(4); hold on;
    set(4,'renderer','zbuffer');
    for ii = 1:3
        plot(vTimes,nanmean(vSDF(proTrials(:,maxInd) & congMat(:,ii),:),1),'color','g','linewidth',1,'linestyle',congStyle{ii});
        plot(vTimes,nanmean(vSDF(antiTrials(:,maxInd) & congMat(:,ii),:),1),'color','r','linewidth',1,'linestyle',congStyle{ii});
    end
    set(gca,'XLim',vXwind);
    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(4,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(4,sprintf('%s/PA-Cong-StimDir.png',saveDir));
        close(4);
    end

    % Figure 5 = pro and anti at different
    % eccentricities
    figure(5);
    for il = 1:length(proLocs)
        sp(il) = subplot(3,3,spInds(il)); hold on;
        for ie = 1:length(uEccs)
            thesePro = proTrials(:,il) & Task.Eccentricity == uEccs(ie);
            theseAnti = antiTrials(:,il) & Task.Eccentricity == uEccs(ie);

            plot(vTimes,nanmean(vSDF(thesePro,:),1),'color',[.1 .9 .1]+[.05 0 .05].*uEccs(ie));
            plot(vTimes,nanmean(vSDF(theseAnti,:),1),'color',[.9 .1 .1]+[0 .05 .05].*uEccs(ie));
        end
        set(gca,'XLim',vXwind);
    end
    linkaxes(sp,'y');
    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(5,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(5,sprintf('%s/PA-Ecc-Stim.png',saveDir));
        close(5);
    end

    % Figure 6 = pro and anti at different
    % eccentricities, saccade aligned
    figure(6);
    for il = 1:length(proLocs)
        sp(il) = subplot(3,3,spInds(il)); hold on;
        for ie = 1:length(uEccs)
            thesePro = proTrials(:,il) & Task.Eccentricity == uEccs(ie);
            theseAnti = antiTrials(:,il) & Task.Eccentricity == uEccs(ie);

            plot(mTimes,nanmean(mSDF(thesePro,:),1),'color',[.1 .9 .1]+[.05 0 .05].*uEccs(ie));
            plot(mTimes,nanmean(mSDF(theseAnti,:),1),'color',[.9 .1 .1]+[0 .05 .05].*uEccs(ie));
        end
        set(gca,'XLim',mXwind);
    end
    linkaxes(sp,'y');
    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(6,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(6,sprintf('%s/PA-Ecc-Sacc.png',saveDir));
        close(6);
    end
    % Here we plot the SDFs the way they are in S&S 03 to demonstrate SST, EST, SRT
    % Get receptive field from pro trials
    for il = 1:length(proLocs)
        vResp(il) = nanmean(nanmean(vSDF(isPro & isCorrect & isPA & Task.TargetLoc==proLocs(il),ismember(vTimes,vWind))));
    end
    [~,maxInd] = max(vResp);
    oppInd = mod(maxInd+length(proLocs)/2,length(proLocs)); if oppInd == 0, oppInd = length(proLocs); end

    % Get the pro/anti specific timing values
    [sstp, proDiff] = getSST(nanmean(vSDF(isCorrect == 1 & isPA & isPro & Task.TargetLoc==proLocs(maxInd),:),1),...
        nanmean(vSDF(isCorrect == 1 & isPA & isPro & Task.TargetLoc==proLocs(oppInd),:),1),vTimes);
    [ssta, antiDiff] = getSST(nanmean(vSDF(isCorrect == 1 & isPA & isAnti & Task.TargetLoc==proLocs(maxInd),:),1),...
        nanmean(vSDF(isCorrect == 1 & isPA & isAnti & Task.TargetLoc==proLocs(oppInd),:),1),vTimes);
    esta = getSST(nanmean(vSDF(isCorrect == 1 & isPA & isAnti & Task.TargetLoc==proLocs(oppInd),:),1),...
        nanmean(vSDF(isCorrect == 1 & isPA & isAnti & Task.TargetLoc==proLocs(maxInd),:),1),vTimes);
    srtp = getSST(proDiff,antiDiff,vTimes);
    srta = getSST(antiDiff,proDiff,vTimes);
    if ~isnan(srtp) && ~isnan(srta)
        srt = min([srtp,srta]);
    elseif isnan(srtp) && ~isnan(srta)
        srt = srta;
    elseif ~isnan(srtp) && isnan(srta)
        srt = srtp;
    else
        srt = nan;
    end

    % Get them plotted
    figure(7);
    % SR = stimulus raw
    sr(1) = subplot(2,3,1); hold on;
    plot(vTimes,nanmean(vSDF(isCorrect == 1 & isPA & isPro & Task.TargetLoc==proLocs(maxInd),:),1),'k','linewidth',3);
    plot(vTimes,nanmean(vSDF(isCorrect == 1 & isPA & isPro & Task.TargetLoc==proLocs(oppInd),:),1),'k','linewidth',1);
    v=vline(sstp); set(v,'color','k','linestyle','--');
    sr(2) = subplot(2,3,4); hold on;
    plot(vTimes,nanmean(vSDF(isCorrect == 1 & isPA & isAnti & Task.TargetLoc==proLocs(maxInd),:),1),'r','linewidth',3);
    plot(vTimes,nanmean(vSDF(isCorrect == 1 & isPA & isAnti & Task.TargetLoc==proLocs(oppInd),:),1),'r','linewidth',1);
    v=vline(ssta); set(v,'color','k','linestyle','--');
    v=vline(esta); set(v,'color','k','linestyle','--');
    linkaxes(sr,'xy');
    yLims = get(sr(1),'YLim');
    axes(sr(1));
    scatter(sstp,yLims(2).*.75,'ok');
    axes(sr(2));
    scatter(ssta,yLims(2).*.75,'ok');
    scatter(esta,yLims(2).*.75,'xk');
    set(sr,'tickdir','out','box','off','XLim',[0,300]);

    % SD = stimulus diff
    diffBlMean = nanmean(proDiff(vTimes >= -50 & vTimes <= 50));
    diffBlStd = nanstd(proDiff(vTimes >= -50 & vTimes <= 50));
    sd(1) = subplot(2,3,2); hold on;
    plot(vTimes,proDiff,'.k');
    blm = hline(diffBlMean); set(blm,'color','k');
    blsp = hline(diffBlMean+2*diffBlStd); set(blsp,'color',[.4 .4 .4],'linestyle','-');
    blsm = hline(diffBlMean-2*diffBlStd); set(blsm,'color',[.4 .4 .4],'linestyle','-');
    blspp = hline(diffBlMean+5*diffBlStd); set(blspp,'color',[.4 .4 .4],'linestyle','--');
    blsmm = hline(diffBlMean-5*diffBlStd); set(blsmm,'color',[.4 .4 .4],'linestyle','--');
    v=vline(sstp); set(v,'color','k','linestyle','--');
    sd(2) = subplot(2,3,5); hold on;
    diffBlMean = nanmean(antiDiff(vTimes >= -50 & vTimes <= 50));
    diffBlStd = nanstd(antiDiff(vTimes >= -50 & vTimes <= 50));
    plot(vTimes,antiDiff,'.r');
    blm = hline(diffBlMean); set(blm,'color','k');
    blsp = hline(diffBlMean+2*diffBlStd); set(blsp,'color',[.4 .4 .4],'linestyle','-');
    blsm = hline(diffBlMean-2*diffBlStd); set(blsm,'color',[.4 .4 .4],'linestyle','-');
    blspp = hline(diffBlMean+5*diffBlStd); set(blspp,'color',[.4 .4 .4],'linestyle','--');
    blsmm = hline(diffBlMean-5*diffBlStd); set(blsmm,'color',[.4 .4 .4],'linestyle','--');
    linkaxes(sd,'xy');
    v=vline(ssta); set(v,'color','k','linestyle','--');
    v=vline(esta); set(v,'color','k','linestyle','--');
    yLims = get(sd(1),'YLim');
    axes(sd(1));
    scatter(sstp,yLims(2).*.75,'ok');
    axes(sd(2));
    scatter(ssta,yLims(2).*.75,'ok');
    scatter(esta,yLims(2).*.75,'xk');
    set(sd,'tickdir','out','box','off','XLim',[0,300]);



    % st = stimulus task
    st(1) = subplot(2,3,3); hold on;
    plot(vTimes,proDiff,'.k');
    plot(vTimes,antiDiff,'.r');
    v=vline(srt); set(v,'color',[.2 .8 .8],'linewidth',3);
    st(2) = subplot(2,3,6); hold on;
    diffTask = proDiff-antiDiff;
    taskBL = nanmean(diffTask(vTimes >= -50 & vTimes <= 50));
    taskStd = nanstd(diffTask(vTimes >= -50 & vTimes <= 50));
    plot(vTimes,diffTask,'-k');
    blm = hline(taskBL); set(blm,'color','k');
    blsp = hline(taskBL+2*taskStd); set(blsp,'color',[.4 .4 .4],'linestyle','-');
    blsm = hline(taskBL-2*taskStd); set(blsm,'color',[.4 .4 .4],'linestyle','-');
    blspp = hline(taskBL+5*taskStd); set(blspp,'color',[.4 .4 .4],'linestyle','--');
    blsmm = hline(taskBL-5*taskStd); set(blsmm,'color',[.4 .4 .4],'linestyle','--');
    v=vline(srt); set(v,'color',[.2 .8 .8],'linewidth',3);
    linkaxes(st,'xy'); set(st,'tickdir','out','box','off','XLim',[0,300]);

    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(7,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(7,sprintf('%s/PA-SST-EST-SRT.png',saveDir));
        close(7);
    end
    
    % Figure 8: Vertical Stimulus in RF as correct pro (red solid),
    % incorrect pro (red dashed), selected incorrect distractor (blue solid),
    % unselected distractor (blue dashed), correct anti distractor (green
    % solid) and incorrect congruent anti (green dashed)
    rfDiff = nan(size(Task.StimDiff,1),1);
    isShapeSing = false(size(Task.StimDiff,1),1);
    for it = 1:size(Task.StimDiff,1)
        if any(~isnan(Task.StimLoc(it,:)))
            rfDiff(it) = Task.StimDiff(it,Task.StimLoc(it,:)==vRF);
            isShapeSing(it) = shapeSingMat(it,Task.StimLoc(it,:)==vRF);
        end
    end
    
    figure(8); hold on;
    plot(vTimes,nanmean(vSDF(isPro & isCorrect & Task.TargetLoc==vRF,:),1),'color',[.8 .2 .2],'linewidth',2);
    plot(vTimes,nanmean(vSDF(isPro & ~isCorrect & isGood & Task.TargetLoc==vRF,:),1),'color',[.8 .2 .2],'linestyle','--','linewidth',2);
    plot(vTimes,nanmean(vSDF(Task.EndStimLoc == vRF & isAnti & isCorrect & rfDiff==0,:),1),'color',[.2 .8 .2],'linewidth',2);
    plot(vTimes,nanmean(vSDF(isAnti & ~isCorrect & isGood & rfDiff == 0 & mod(Task.TargetLoc+180,360)==vRF,:),1),'color',[.2 .8 .2],'linestyle','--','linewidth',2);
    plot(vTimes,nanmean(vSDF(Task.TargetLoc ~= vRF & rfDiff==0 & Task.EndStimLoc==vRF,:),1),'color',[.2 .2 .8],'linewidth',2);
    plot(vTimes,nanmean(vSDF(Task.TargetLoc ~= vRF & rfDiff==0 & Task.EndStimLoc~=vRF & isGood,:),1),'color',[.2 .2 .8],'linestyle','--','linewidth',2);
    set(gca,'XLim',[-100,300]);
    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(8,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(8,sprintf('%s/PA-VertInRF-Vis.png',saveDir));
        close(8);
    end
    
    figure(9); hold on;
    plot(mTimes,nanmean(mSDF(isPro & isCorrect & Task.TargetLoc==vRF,:),1),'color',[.8 .2 .2],'linewidth',2);
    plot(mTimes,nanmean(mSDF(isPro & ~isCorrect & isGood & Task.TargetLoc==vRF,:),1),'color',[.8 .2 .2],'linestyle','--','linewidth',2);
    plot(mTimes,nanmean(mSDF(Task.EndStimLoc == vRF & isAnti & isCorrect & rfDiff==0,:),1),'color',[.2 .8 .2],'linewidth',2);
    plot(mTimes,nanmean(mSDF(isAnti & ~isCorrect & isGood & rfDiff == 0 & mod(Task.TargetLoc+180,360)==vRF,:),1),'color',[.2 .8 .2],'linestyle','--','linewidth',2);
    plot(mTimes,nanmean(mSDF(Task.TargetLoc ~= vRF & rfDiff==0 & Task.EndStimLoc==vRF,:),1),'color',[.2 .2 .8],'linewidth',2);
    plot(mTimes,nanmean(mSDF(Task.TargetLoc ~= vRF & rfDiff==0 & Task.EndStimLoc~=vRF & isGood,:),1),'color',[.2 .2 .8],'linestyle','--','linewidth',2);
    set(gca,'XLim',[-300,100]);
    suptitle(sprintf('%s-%s',sessName,analyzeChannels{ic}));
    if doSave
        set(9,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(9,sprintf('%s/PA-VertInRF-Mov.png',saveDir));
        close(9);
    end
    
    
    % Figure 10: V vs H vs S, not selected, in RF
    figure(10); 
    sp(1) = subplot(1,2,1); hold on;
    plot(vTimes,nanmean(vSDF(rfDiff==0 & Task.EndStimLoc ~= vRF & Task.TargetLoc ~= vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.8 .2 .2]);
    plot(vTimes,nanmean(vSDF(rfDiff==1 & Task.EndStimLoc ~= vRF & Task.TargetLoc ~= vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.2 .8 .2]);
    plot(vTimes,nanmean(vSDF(rfDiff==2 & Task.EndStimLoc ~= vRF & Task.TargetLoc ~= vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.2 .2 .8]);
    set(gca,'XLim',[-200,300]);
    sp(2) = subplot(1,2,2); hold on;
    plot(mTimes,nanmean(mSDF(rfDiff==0 & Task.EndStimLoc ~= vRF & Task.TargetLoc ~= vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.8 .2 .2]);
    plot(mTimes,nanmean(mSDF(rfDiff==1 & Task.EndStimLoc ~= vRF & Task.TargetLoc ~= vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.2 .8 .2]);
    plot(mTimes,nanmean(mSDF(rfDiff==2 & Task.EndStimLoc ~= vRF & Task.TargetLoc ~= vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.2 .2 .8]);
    set(gca,'XLim',[-300,200]);
    linkaxes(sp,'y');
    
    % Figure 11: V in RF selected vs V in RF not selected (off distractors)
    figure(11);
    sp(1) = subplot(1,2,1); hold on;
    plot(vTimes,nanmean(vSDF(rfDiff==0 & Task.EndStimLoc==vRF & Task.TargetLoc~=vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.2 .8 .2]);
    plot(vTimes,nanmean(vSDF(rfDiff==0 & Task.EndStimLoc~=vRF & Task.TargetLoc~=vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.8 .2 .2]);
    set(gca,'XLim',[-200,300]);
    sp(1) = subplot(1,2,2); hold on;
    plot(mTimes,nanmean(mSDF(rfDiff==0 & Task.EndStimLoc==vRF & Task.TargetLoc~=vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.2 .8 .2]);
    plot(mTimes,nanmean(mSDF(rfDiff==0 & Task.EndStimLoc~=vRF & Task.TargetLoc~=vRF & mod(Task.TargetLoc+180,360)~=vRF,:),1),'color',[.8 .2 .2]);
    set(gca,'XLim',[-300,200]);
    linkaxes(sp,'y');
    
    % Figure 12: V Correctly Selected vs V Incorrectly Selected
    figure(12);
    sp(1) = subplot(1,2,1); hold on;
    plot(vTimes,nanmean(vSDF(rfDiff==0 & Task.EndStimLoc==vRF & isCorrect & isPro,:),1),'color','k');
    plot(vTimes,nanmean(vSDF(rfDiff==0 & Task.EndStimLoc==vRF & ~isCorrect & isGood & isPro,:),1),'color','k','linestyle','--');
    plot(vTimes,nanmean(vSDF(rfDiff==0 & Task.EndStimLoc==vRF & isCorrect & isAnti,:),1),'color',[.8 .2 .2]);
    plot(vTimes,nanmean(vSDF(rfDiff==0 & Task.EndStimLoc==vRF & ~isCorrect & isGood & isAnti,:),1),'color',[.8 .2 .2],'linestyle','--');
    set(gca,'XLim',[-200,300]);
    sp(2) = subplot(1,2,2); hold on;
    plot(mTimes,nanmean(mSDF(rfDiff==0 & Task.EndStimLoc==vRF & isCorrect & isPro,:),1),'color','k');
    plot(mTimes,nanmean(mSDF(rfDiff==0 & Task.EndStimLoc==vRF & ~isCorrect & isGood & isPro,:),1),'color','k','linestyle','--');
    plot(mTimes,nanmean(mSDF(rfDiff==0 & Task.EndStimLoc==vRF & isCorrect & isAnti,:),1),'color',[.8 .2 .2]);
    plot(mTimes,nanmean(mSDF(rfDiff==0 & Task.EndStimLoc==vRF & ~isCorrect & isGood & isAnti,:),1),'color',[.8 .2 .2],'linestyle','--');
    set(gca,'XLim',[-300,200]);
    linkaxes(sp,'y');
    
    % Figure 13: Color Singleton Only vs Shape Singleton Only vs Both vs
    % Neither
    figure(13);
    sp(1) = subplot(1,2,1); hold on;
    % Color Singleton
    plot(vTimes,nanmean(vSDF(Task.TargetLoc==vRF & ~isShapeSing & Task.EndStimLoc==vRF,:),1),'color',[.8 .2 .2]);
    % Color + Shape
    plot(vTimes,nanmean(vSDF(Task.TargetLoc==vRF & isShapeSing & Task.EndStimLoc==vRF,:),1),'color',[.8 .2 .2],'linestyle','--');
    % Shape only
    plot(vTimes,nanmean(vSDF(Task.TargetLoc~=vRF & isShapeSing & Task.EndStimLoc==vRF,:),1),'color','k','linestyle','--');
    % Neither
    plot(vTimes,nanmean(vSDF(Task.TargetLoc~=vRF & ~isShapeSing & Task.EndStimLoc==vRF,:),1),'color','k');
    set(gca,'XLim',[-200,300]);
    sp(2) = subplot(1,2,2); hold on;
    % Color Singleton
    plot(mTimes,nanmean(mSDF(Task.TargetLoc==vRF & ~isShapeSing & Task.EndStimLoc==vRF,:),1),'color',[.8 .2 .2]);
    % Color + Shape
    plot(mTimes,nanmean(mSDF(Task.TargetLoc==vRF & isShapeSing & Task.EndStimLoc==vRF,:),1),'color',[.8 .2 .2],'linestyle','--');
    % Shape only
    plot(mTimes,nanmean(mSDF(Task.TargetLoc~=vRF & isShapeSing & Task.EndStimLoc==vRF,:),1),'color','k','linestyle','--');
    % Neither
    plot(mTimes,nanmean(mSDF(Task.TargetLoc~=vRF & ~isShapeSing & Task.EndStimLoc==vRF,:),1),'color','k');
    set(gca,'XLim',[-300,200]);
    linkaxes(sp,'y');
    
    
    if pause
        keyboard
    end
    
    

end

end


    
function [sst, diffSDF] = getSST(in,across,times)
    % Get difference and z-score
    diffSDF = in - across;
    blMean = nanmean(diffSDF(times >= -50 & times <= 50));
    blStd = nanstd(diffSDF(times >= -50 & times <= 50));
    zSDF = (diffSDF-blMean)./blStd;
    
    % Get the first time that zSDF >= 2 for more than 15ms
    overFive = klGetConsecutive(zSDF >= 5);
    overTwo = klGetConsecutive(zSDF >= 2);
    firstOver = find(overFive >= 15 & times >= 0,1);
    if overTwo(firstOver) >= 15
        sst = times(find(overTwo(1:firstOver) == 0,1,'last'));
    else
        sst = nan;
    end    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    