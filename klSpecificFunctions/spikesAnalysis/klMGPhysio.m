function returned = klMGPhysio(sessName,varargin)

baseDir = [tebaMount,'Users/Kaleb/proAntiProcessed/';
doSave = 1;
watch = 0;
fresh = 1;
channel = 0;
visWind = 50:150;
movWind = -50:0;
fivePos = [.4108, .4096 .2027 .2157];

varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-s'}
            doSave = varargin{varStrInd(iv)+1};
        case {'-f'}
            fresh = varargin{varStrInd(iv)+1};
        case {'-c'}
            channel = varargin{varStrInd(iv)+1};
    end
end

% Load in behavior and eye positions
% sessDir = chanStruct.folder;
% [~,sessName] = fileparts(sessDir);
sessDir = [baseDir,filesep,sessName];
% [~,sessName] = fileparts(sessDir);
load([sessDir, '/Behav.mat']);

% Find unit .mat files
% chanMats = dir([sessDir,'/chan*.mat']);
allUnits = dir([sessDir,filesep,'Chan*']);
unitNames = {allUnits.name};
if channel
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
colors = 'rgbcmk';

% Make some logicals
isMG = strcmpi(Task.TaskType,'MG');
isCorrect = Task.Correct == 1;
isGood = Task.Abort == 0;

% get eccentricities and target locations
uEccs = nunique(Task.Eccentricity(isMG));
uLocs = 0:45:315;
polLocs = klDeg2Rad(uLocs); polLocs(end+1) = polLocs(1);
spInds = [6 3 2 1 4 7 8 9];


keepVars = whos; keepVars = {keepVars.name}; keepVars{end+1} = 'keepVars'; keepVars{end+1} = 'clrVars'; keepVars{end+1} = 'ic';

% Start unit loop
for ic = 1:length(analyzeChannels)
    close all;
    currVars = whos; currVars = {currVars.name};
    clrVars = currVars(~ismember(currVars,keepVars));
    for iv = 1:length(clrVars), clear(clrVars{iv}); end
    
    fprintf('Working on unit %d (of %d)...\n',ic,length(analyzeChannels));
    saveDir = [sessDir,filesep,analyzeChannels{ic}];
    if doSave && length(dir(sprintf('%s/MG*.png',saveDir))) == 2 && ~fresh
        continue
    end
    
    
    % load([chanStruct.folder,filesep,chanStruct.name]);
    load([sessDir,filesep,analyzeChannels{ic},filesep,lower(analyzeChannels{ic}(1)),analyzeChannels{ic}(2:end),'.mat']);
    spkMat = klPlaceEvents(Task,spkTimes);
    
    if size(spkMat,2) < 5
        returned = 1;
        continue
    end
    returned = 0;
    
    [vSDF,vTimes] = klSpkRatev2(spkMat,'-q',1);
    [mSDF,mTimes] = klSpkRatev2(spkMat-repmat(Task.SRT+Task.GoCue,1,size(spkMat,2)),'-q',1);
    
    
    % figure 1: Vis responses
    figure(1);
    for il = 1:length(uLocs)
        sp(il) = subplot(3,3,spInds(il));
        for ie = 1:length(uEccs)
            myTrials = isMG & isCorrect & Task.Eccentricity == uEccs(ie) & Task.TargetLoc == uLocs(il);
            if sum(myTrials) >= 5
                plot(vTimes,nanmean(vSDF(myTrials,:),1),'color',colors(ie)); hold on;
            end
            mnVisActivity(ie,il) = nanmean(nanmean(vSDF(myTrials,ismember(vTimes,visWind))));
        end
        set(gca,'XLim',vXwind);
    end
    linkaxes(sp,'y');
%     pp=subplot(3,3,5,polaraxes);
    pp = polaraxes('Position',fivePos);
    mnVisActivity(:,end+1) = mnVisActivity(:,1);
    for ie = 1:length(uEccs)
        axes(pp);
        hp(ie) = polarplot(polLocs,mnVisActivity(ie,:)); hold on;
        set(hp(ie),'color',colors(ie),'linewidth',1.25);
    end
    set(pp,'ThetaTickLabel',[],'RTickLabel',[]);
    suptitle('Vis Responses');
    if doSave
        set(1,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(1,sprintf('%s/MG-VisResponses.png',saveDir));
        close(1);
    end
    clear sp hp pp
    pause(.2);
    
    % figure 2: mov responses
    figure(2);
    for il = 1:length(uLocs)
        sp(il) = subplot(3,3,spInds(il));
        for ie = 1:length(uEccs)
            myTrials = isMG & isCorrect & Task.Eccentricity == uEccs(ie) & Task.TargetLoc == uLocs(il);
            if sum(myTrials) >= 5
                plot(mTimes,nanmean(mSDF(myTrials,:),1),'color',colors(ie)); hold on;
            end
            mnMovActivity(ie,il) = nanmean(nanmean(mSDF(myTrials,ismember(mTimes,movWind))));
        end
        set(gca,'XLim',mXwind);
    end
    linkaxes(sp,'y');
    mnMovActivity(:,end+1) = mnMovActivity(:,1);
    pp = polaraxes('Position',fivePos);
    for ie = 1:length(uEccs)
        axes(pp);
        hp(ie)=polarplot(polLocs,mnMovActivity(ie,:)); hold on;
        set(hp(ie),'color',colors(ie),'linewidth',1.25);
    end
    set(pp,'ThetaTickLabel',[],'RTickLabel',[]);
    suptitle('Mov Responses');
    if doSave
        set(2,'paperposition',[.2 .1 10.5 7.5],'papersize',[11,8]);
        saveas(2,sprintf('%s/MG-MovResponses.png',saveDir));
        close(2);
    end
    
end























