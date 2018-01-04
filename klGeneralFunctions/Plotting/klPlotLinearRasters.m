function h = klPlotLinearRasters(sessFold,recSys,varargin)


% Set defaults
chanSkip = 200;

% Do the silly global thing
sessOld = sessFold;
recOld = recSys;

global tebaDir sessFold recSys Task monk

switch recSys
    case {'tdt'}
        tebaDir = '/mnt/teba/Users/Kaleb/dataProcessed';
    case {'plex'}
        tebaDir = '/mnt/teba/Users/Wolf/ephys_db';
end

sessFold = sessOld;
recSys = recOld;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-t'}
            tebaDir = varargin{varStrInd(iv)+1};
        case {'-m'}
            monk = varargin{varStrInd(iv)+1};
    end
end



% Let's start by ending this function and making a data adapter below

% Get matFile struct
chanStruct = dataAdapt;

% Loop through and concatenate x and spikes
vSpikes = [];
mSpikes = [];
vY = [];
mY = [];
offset = 0;
for ic = 1:length(chanStruct)
    % Start from the bottom...
    vSpikes = cat(2,vSpikes,chanStruct(ic).vSpks);
    vY = cat(2,vY,repmat((1:size(chanStruct(ic).vSpks,1))'+offset,1,size(chanStruct(ic).vSpks,2)));
    mSpikes = cat(2,mSpikes,chanStruct(ic).mSpks);
    mY = cat(2,mY,repmat((1:size(chanStruct(ic).mSpks,1))'+offset,1,size(chanStruct(ic).mSpks,2)));
    if offset & mod(ic,32)==1
        offset = 0;
        % Open the figures
        vFig = figure();
        hv = scatter(vSpikes(:),vY(:),[],'k','marker','.');
        set(gca,'XLim',[-100,250]);
        mFig = figure();
        hm = scatter(mSpikes(:),mY(:),[],'k','marker','.');
        set(gca,'XLim',[-250,100]);
        vSpikes = [];
        mSpikes = [];
        vY = [];
        mY = [];

    else        
        offset = offset + size(chanStruct(ic).mSpks,1)+chanSkip;
    end
    
end

  
    

% Open the figures
% vFig = figure();
% sp(1) = subplot(1,2,1);
% hv = scatter(vSpikes(:),vY(:),[],[.4 .4 .4],'marker','.');
% set(gca,'box','off','XLim',[-50,300],'tickdir','out','ticklength',get(gca,'ticklength').*4);
% xlabel('Time From Stimulus'); ylabel('Trials by Channel');
% sp(2) = subplot(1,2,2);
% hm = scatter(mSpikes(:),mY(:),[],[.4 .4 .4],'marker','.');
% set(gca,'box','off','XLim',[-150,100],'tickdir','out','ticklength',get(gca,'ticklength').*4,'yaxisloc','right');
% xlabel('Time From Saccade'); ylabel('Trials by Channel');
% linkaxes(sp,'y');
% set(sp,'YLim',[0,offset]);
% suptitle([sessFold,' Rasters (A3L2)']);
% keyboard

figure();
nOff = 0;
sp(1) = subplot(1,2,1);
sp(2) = subplot(1,2,2);
for ic = length(chanStruct):-1:1
    axes(sp(1));
    if ~isempty(chanStruct(ic).vSDF)
        vScale = chanStruct(ic).vSDF./max(abs(nanmean(chanStruct(ic).vSDF(ismember(chanStruct(ic).vTimes,[-50:250])))));
        bl = nanmean(nanmean(vScale(:,ismember(chanStruct(ic).vTimes,-300:-100)),1),2);
        plot(chanStruct(ic).vTimes,nanmean(vScale,1)+nOff-bl,'k'); hold on;
    end
    axes(sp(2));
    if ~isempty(chanStruct(ic).mSDF)
        mScale = chanStruct(ic).mSDF./max(abs(nanmean(chanStruct(ic).mSDF(ismember(chanStruct(ic).mTimes,[-200:100])))));
        plot(chanStruct(ic).mTimes,nanmean(mScale,1)+nOff-bl,'k'); hold on;
    end
    nOff = nOff+1;
end
set(sp(1),'XLim',[-50,250],'box','off','tickdir','out','ticklength',get(sp(1),'ticklength').*3);
set(sp(2),'XLim',[-200,100],'box','off','yaxisloc','right','tickdir','out','ticklength',get(sp(2),'ticklength').*3);
suptitle(sessFold);
% Start with vis
% figure(vFig);


end


function chanStruct = dataAdapt()

global sessFold recSys Task tebaDir monk

switch recSys
    case {'tdt'}
        % Load in Behav.mat
        load([tebaDir,filesep,sessFold,filesep,'Behav.mat']);
        
        % Get channel .mats
        chanFolds = dir([tebaDir,filesep,sessFold,filesep,'Channel*']);
        chanNames = {chanFolds.name};
        hasMUA = cellfun(@(x) exist([sessFold,filesep,x,filesep,'chan*MUA.mat']),chanNames);
        whichChan = cellfun(@(x) str2double(x(8:end)),chanNames);
        for ic = 1:length(whichChan)
            if exist('printStr','var')
                fprintf(repmat('\b',1,length(printStr)));
            end
            printStr = sprintf('Getting Channel %d...\n',ic);
            fprintf(printStr);
            chan = whichChan(ic);
            try
                load([tebaDir,filesep,sessFold,filesep,'Channel',num2str(chan),filesep,'chan',num2str(chan),'MUA.mat']);
                spkMat = klPlaceEvents(Task,spkTimes);
                isMG = strcmpi(Task.TaskType,'MG');
                spkMat = spkMat(isMG,:);
                taskMG = klCutTask(Task,isMG);
                cutAgain = taskMG.Correct==1 & ismember(Task.TargetLoc,[315,0,45]);
                taskMG = klCutTask(taskMG,cutAgain);
                spkMat = spkMat(cutAgain,:);
                [mSDF,mTimes] = klSpkRatev2(spkMat-taskMG.SRT+taskMG.GoCue,'-q',1);
                mSpks = spkMat-taskMG.SRT+taskMG.GoCue;
                cSpks = spkMat;
                cSpks(cSpks >= repmat(taskMG.SRT+taskMG.GoCue,1,size(cSpks,2))) = nan;
                [vSDF,vTimes] = klSpkRatev2(cSpks,'-q',1);
                chanStruct(chan).vSDF = vSDF;
                chanStruct(chan).vSpks = cSpks;
                chanStruct(chan).vTimes = vTimes;
                chanStruct(chan).mSDF = mSDF;
                chanStruct(chan).mSpks = mSpks;
                chanstruct(chan).mTimes = mTimes;
            catch
                chanStruct(chan).vSDF = [];
                chanStruct(chan).vSpks = [];
                chanStruct(chan).vTimes = [];
                chanStruct(chan).mSDF = [];
                chanStruct(chan).mSpks = [];
                chanStruct(chan).mTimes = [];
            end
        end
    case {'plex'}
        % Get channel folders
        chanFolds = dir([tebaDir,filesep,monk,filesep,sessFold,filesep,'DSP/DSP*']);
        chanNames = {chanFolds.name};
        chanNums = cellfun(@(x) str2double(x(4:5)),chanNames);
        
        for ic = 1:max(chanNums)
            chanStruct(ic).vSpks = [];
            chanStruct(ic).mSpks = [];
            fprintf('Getting Channel %d...\n',ic);
            unitThisChan = find(chanNums==ic);
            if ~isempty(unitThisChan)
                for ii = 1:length(unitThisChan)
                    mgFile = dir([tebaDir,filesep,monk,filesep,sessFold,filesep,'DSP',filesep,chanNames{unitThisChan(ii)},'/*MG*.mat']);
                    if ~isempty(mgFile)
                        load([mgFile(1).folder,filesep,mgFile(1).name]);
                        chanStruct(ic).mSpks = cat(2,chanStruct(ic).mSpks,spiketimes-repmat(Task.SRT+Task.GoCue,1,size(spiketimes,2)));
                        spiketimes(spiketimes >= repmat(Task.SRT+Task.GoCue,1,size(spiketimes,2))) = nan;
                        chanStruct(ic).vSpks = cat(2,chanStruct(ic).vSpks,spiketimes);
                    end
                end
                chanStruct(ic).vSpks = sort(chanStruct(ic).vSpks,2);
                chanStruct(ic).mSpks = sort(chanStruct(ic).mSpks,2);
                [chanStruct(ic).vSDF,chanStruct(ic).vTimes] = klSpkRatev2(chanStruct(ic).vSpks);
                [chanStruct(ic).mSDF,chanStruct(ic).mTimes] = klSpkRatev2(chanStruct(ic).mSpks);
            end
        end
end
                
                
            
        
        
        
        
        
end

            
            
         