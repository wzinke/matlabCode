function returned = klMGPhysio(Task,spkTimes,varargin)

baseDir = [tebaMount,'Users/Kaleb/proAntiProcessed/';
doSave = 1;
watch = 0;
fresh = 1;

varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-s'}
            doSave = varargin{varStrInd(iv)+1};
        case {'-f'}
            fresh = varargin{varStrInd(iv)+1};
        case {'-c'}
            chanBase = varargin{varStrInd(iv)+1};
    end
end

if ~exist('chanBase','var')
    chanBase = 1;
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
spInds = [6 3 2 1 4 7 8 9];

% load([chanStruct.folder,filesep,chanStruct.name]);
spkMat = klPlaceEvents(Task,spkTimes);

if size(spkMat,2) < 5
    returned = 1;
    return
end
returned = 0;

[vSDF,vTimes] = klSpkRatev2(spkMat,'-q',1);
[mSDF,mTimes] = klSpkRatev2(spkMat-repmat(Task.SRT+Task.GoCue,1,size(spkMat,2)),'-q',1);


% figure 1: Vis responses
figure(chanBase); set(gcf,'Position',[562.3333  765.0000  560.0000  420.0000]);
for il = 1:length(uLocs)
    sp(il) = subplot(3,3,spInds(il));
    for ie = 1:length(uEccs)
        myTrials = isMG & isGood & Task.Eccentricity == uEccs(ie) & Task.TargetLoc == uLocs(il);
        if sum(myTrials) >= 5
            plot(vTimes,nanmean(vSDF(myTrials,:),1),'color',colors(ie)); hold on;%,'linewidth',ie); hold on;
        end
    end
    set(gca,'XLim',vXwind);
end
linkaxes(sp,'y');
suptitle('Vis Responses');

% figure 2: mov responses
figure(chanBase + 100); set(gcf,'Position',[1.1703    0.7630    0.5600    0.4200].*1e3);
for il = 1:length(uLocs)
    sp(il) = subplot(3,3,spInds(il));
    for ie = 1:length(uEccs)
        myTrials = isMG & isCorrect & Task.Eccentricity == uEccs(ie) & Task.TargetLoc == uLocs(il);
        if sum(myTrials) >= 5
            plot(mTimes,nanmean(mSDF(myTrials,:),1),'color',colors(ie)); hold on; %,'linewidth',ie); hold on;
        end
    end
    set(gca,'XLim',mXwind);
end
linkaxes(sp,'y');
suptitle('Mov Responses');



    






















