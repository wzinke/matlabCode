function [allSDF,allSDFTimes,allAreas,allRows,allMonks,allOtherTasks,allSess,allChans,allParams] = klPullPlexSDFsv2(task,varargin)


switch task
    case {'MG','mg'}
        xlFile = 'klDataBookKeeping_mg.xlsx';
        headRow = 4;
        sessCol = 2;
        chanCol = 5;
    case {'Search','Capture','cap','search','Cap'}
        xlFile = 'cosDataBookKeeping_capturev2.xlsx';
        headRow = 3;
        sessCol = 1;
        chanCol = 2;
end

% Set data location
baseDir = [tebaMount,'Users/Wolf/ephys_db/'];
xlFile = [mlRoot,'/Dropbox/Schall-Lab/dataExcels/plexBookKeeping/',xlFile];

% Set constants
minRate = 0;
isiCrit = inf;
snrCrit = .85;
areaCrit = {'FEF','F2'};
doReward = 0;
kType = 'psp';
gWidth = 20;
clip = 0;
waveTimes = ((1:32)-9).*25;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-a','area'}
            areaCrit = varargin{varStrInd(iv)+1};
        case {'-r','reward'}
            doReward = varargin{varStrInd(iv)+1};
        case {'-k'}
            kType = varargin{varStrInd(iv)+1};
        case {'-w'}
            gWidth = varargin{varStrInd(iv)+1};
        case {'-c'}
            clip = varargin{varStrInd(iv)+1};
    end
end

% Set helpful windows
vWind = -300:500;
mWind = -500:300;
rWind = -200:200;
blWind = -300:-200;
vCheck = 50:150;
mCheck = -50:0;

if doReward
    nAligns = 3;
else
    nAligns = 2;
end

monks = {'Gauss','Helmholtz','Darwin'};
allSDF = {};
allSDFTimes = {};
allAreas = {};
allMonks = {};
allRows = [];
allOtherTasks = {};
allSess = {};
allChans = {};
for im = 1:length(monks)
    fprintf('Pulling data from monkey %s...\n',monks{im});
    monkFolder = [baseDir,monks{im}];
    
    % Read in excel file
    [num,text,all]=xlsread(xlFile,monks{im});
    
    % Get BL Firing Rate column, snr column, and isi column
    snrCol=find(cellfun(@(x) strcmpi(x,'SNR'),all(headRow,:)),1);
    mnRateCol = find(cellfun(@(x) strcmpi(x,'meanRate'),all(headRow,:)),1);
    isiColRate = find(cellfun(@(x) strcmpi(x,' < 3ms'),all(headRow,:)),1);
    areaCol = find(cellfun(@(x) strcmpi(x,'area'),all(headRow,:)),1);
    
    % Use criteria to pick out good rows
    goodRows = cell2mat(all((headRow+1):end,snrCol)) > snrCrit & cell2mat(all((headRow+1):end,mnRateCol)) > minRate & cell2mat(all((headRow+1):end,isiColRate)) < isiCrit & ismember(all((headRow+1):end,areaCol),areaCrit);
    rowInds = find(goodRows)+headRow;
    
    monkSDF = cell(length(rowInds),nAligns);
    monkSDFTimes = cell(length(rowInds),nAligns);
    monkAreas = all(rowInds,areaCol);
%     monkSess = all(rowInds,sessCol)
    monkMonks = cell(length(rowInds),1); [monkMonks{1:length(rowInds)}] = deal(monks{im});
    monkOtherTasks = cell(length(rowInds),1);
    monkSess = all(rowInds,sessCol);
    monkChans = all(rowInds,chanCol);
    
    [monkParams(1:length(rowInds))] = deal(struct('theta',[],'sigma',[],'alpha',[],'offset',[],'blRate',[],'Fano',[],'CV',[],'CV2',[],'LV',[],'LVR',[],'vGOF',[],'mGOF',[],'spkWidth',[]));
    % Start row loop
    for ir = 1:length(rowInds)
        myRow = rowInds(ir);
        printHead = sprintf('Analyzing row %d (%d of %d):  ',myRow,ir,length(rowInds));
        fprintf(printHead);

        % Find and load the appropriate .mat files for spikes and behavior
        printStr = 'Loading data...';
        fprintf(printStr);
        fileBase = sprintf('%s/%s/DSP/%s',monkFolder,all{myRow,sessCol},all{myRow,chanCol});
        switch task
            case {'MG','mg'}
                fileName = all{myRow,4};
                waveFiles = dir([fileBase,'/waves/*MG*.mat']);
                otherCut = 'MG';
            case {'Search','search'}
                fileTmp = dir(sprintf('%s/%s/DSP/%s/*Search*.mat',monkFolder,all{myRow,sessCol},all{myRow,chanCol}));
                if isempty(fileTmp)
                    fileTmp = dir(sprintf('%s/%s/DSP/%s/*Cap*.mat',monkFolder,all{myRow,sessCol},all{myRow,chanCol}));
                end
                if isempty(fileTmp)
                    continue
                else
                    fileName = fileTmp.name;
                end
                waveTmp = dir([fileBase,'/waves/*Search*.mat']);
                if isempty(waveTmp)
                    waveTmp = dir([fileBase,'/waves/*Cap*.mat']);
                end
                waveFiles = waveTmp;
                otherCut = {'Search','Cap','Capture'};
            case {'Capture','cap','Cap','capture'}
                fileTmp = dir(sprintf('%s/%s/DSP/%s/*Cap*.mat',monkFolder,all{myRow,sessCol},all{myRow,chanCol}));
                if isempty(fileTmp)
                    continue
                else
                    fileName = fileTmp.name;
                end
                waveFiles = dir([fileBase,'/waves/*Cap*.mat']);
                otherCut = {'Capture','Cap'};
        end
        fileStr = [fileBase,filesep,fileName];
        
        % Get all tasks here
        cutStr = sprintf('%s_%s_',all{myRow,sessCol},all{myRow,chanCol});
        allTasks = dir(sprintf('%s/%s/DSP/%s/%s*.mat',monkFolder,all{myRow,sessCol},all{myRow,chanCol},cutStr));
        cutStarts = cellfun(@(x) strrep(x,cutStr,''),{allTasks.name},'UniformOutput',0);
        taskStrs = cellfun(@(x) x(1:(strfind(x,'_')-1)),cutStarts,'UniformOutput',0);
        monkOtherTasks{ir} = taskStrs(~ismember(taskStrs,otherCut));
        load(fileStr);
        if ~isempty(waveFiles)
            load(fullfile(waveFiles(1).folder,waveFiles(1).name));
%             wvWidth = klWvWidthv3(spline(1:32,nanmean(wave.waves,1),1:.1:32),spline(1:32,waveTimes,1:.1:32));
            splWave = spline(1:32,nanmean(wave.waves,1),1:.1:32);
            splT = spline(1:32,waveTimes,1:.1:32);
            wvWidth = abs(splT(splWave==max(splWave))-splT(splWave==min(splWave)));
        else
            wvWidth = nan;
        end
        fprintf(repmat('\b',1,length(printStr)));

        printStr = 'Calculating SDFs...';
        fprintf(printStr);
        
        % Calculate SDF for MG and in all MG locations
        spikeMat = spiketimes;%klPlaceEvents(Task,spkTimes);

        % Chop down Task a little bit...
        mgTrials = ismember(Task.TaskType,'MG');
%         Task = klCutTask(Task,mgTrials);
%         spikeMat = spikeMat(mgTrials,:);
        switch task
            case {'MG','mg'}
                goodTrials = ismember(Task.TaskType,'MG');
            case {'Search','search'}
                goodTrials = ismember(Task.TaskType,{'Search','Cap'});
                mgTrials = mgTrials & ismember(Task.TargetLoc,0:90:359);
            case {'Capture','capture','Cap','cap'}
                goodTrials = ismember(Task.TaskType,{'Cap'});
                mgTrials = mgTrials & ismember(Task.TargetLoc,0:90:359);
        end
        spikeMat = spikeMat(goodTrials,:);
        Task = klCutTask(Task,goodTrials);
        
        % Get RF
        rf = getRFv2(spikeMat,Task,'-s',1);
        
        % Get SDFs
        badV = 0; badM = 0; badR = 0;
        mMat = spikeMat-repmat(Task.SRT + Task.GoCue,1,size(spikeMat,2));
        if sum(isfinite(mMat(:))) < 5
            mSDF = 0; mTimes = 0; badM = 1;
        else
            [mSDF,mTimes] = klSpkRatev2(spikeMat-repmat(Task.SRT + Task.GoCue,1,size(spikeMat,2)),'-q',1,'-k',kType,'-w',gWidth);
        end
        rMat = spikeMat-repmat(Task.Reward,1,size(spikeMat,2));
        if sum(isfinite(rMat(:))) < 5
            rSDF = 0; rTimes = 0; badR = 1;
        else
            [rSDF,rTimes] = klSpkRatev2(rMat,'-q',1,'-k',kType,'-w',gWidth);
        end
        if clip
            spikeMat(spikeMat > repmat(Task.SRT+Task.GoCue,1,size(spikeMat,2))) = nan;
        end
        if sum(isfinite(spikeMat(:))) < 5
            vSDF = 0; vTimes = 0; badV = 1;
        else
            [vSDF,vTimes] = klSpkRatev2(spikeMat,'-q',1,'-k',kType,'-w',gWidth);
        end
        
        % Get spike params
        clear vFit mFit
        if badV
            [vFit.mu,vFit.sig,vFit.amp,vFit.bl] = deal(nan);
            vFit.gof = struct('sse',nan,'rsquare',nan,'dfe',nan,'adjrsquare',nan,'rmse',nan);
        else
            vFit = fitGauss(Task,vSDF,vTimes,vCheck);
        end
        if badM
            [mFit.mu,mFit.sig,mFit.amp,mFit.bl] = deal(nan);
            mFit.gof = struct('sse',nan,'rsquare',nan,'dfe',nan,'adjrsquare',nan,'rmse',nan);
        else
            mFit = fitGauss(Task,mSDF,mTimes,mCheck);
        end
        monkParams(ir).theta = [vFit.mu,mFit.mu];
        monkParams(ir).sigma = [vFit.sig,mFit.sig];
        monkParams(ir).alpha = [vFit.amp,mFit.amp];
        monkParams(ir).offset = [vFit.bl,mFit.bl];
        monkParams(ir).vGOF = vFit.gof;
        monkParams(ir).mGOF = mFit.gof;
        monkParams(ir).blRate = nanmean(nanmean(vSDF(:,ismember(vTimes,blWind)),2),1);
        monkParams(ir).CV = klGetCV(spikeMat);
        monkParams(ir).CV2 = klGetCV(spikeMat,'-type','local');
        monkParams(ir).LV = klGetLV(spikeMat);
        monkParams(ir).LVR = klGetLV(spikeMat,'-type','revised');
        monkParams(ir).Fano = klGetFano(spikeMat);
        monkParams(ir).spkWidth = wvWidth;%klWvWidthv3();
        
        fprintf(repmat('\b',1,length(printStr)));
        printStr = 'Putting into output cell...';
        fprintf(printStr);
        % Pull mov and raw visual SDF into monkSDF and monkSDFTimes
        switch task
            case {'MG','mg'}
                if ~badM
                    monkSDF{ir,2} = nanmean(mSDF(Task.TargetLoc==rf(2),ismember(mTimes,mWind)),1);
                    monkSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
                else
                    monkSDF{ir,2} = nan(1,length(mWind));
                    monkSDFTimes{ir,2} = mWind;
                end
%                 if ~badV
%                     monkSDF{ir,3} = nanmean(vSDF(Task.TargetLoc==rf(1),ismember(vTimes,vWind)),1);
%                     monkSDFTimes{ir,3} = vTimes(ismember(vTimes,vWind));
%                 else
%                     monkSDF{ir,3} = nan(1,length(vWind));
%                     monkSDFTimes{ir,3} = vWind;
%                 end
                if doReward
                    if ~badR
                        monkSDF{ir,3} = nanmean(rSDF(ismember(Task.TargetLoc,rf),ismember(rTimes,rWind)),1);
                        monkSDFTimes{ir,3} = rTimes(ismember(rTimes,rWind));
                    else
                        monkSDF{ir,3} = nan(1,length(rWind));
                        monkSDFTimes{ir,3} = rWind;
                    end
                end
                
                % Clip post-mov spikes from vis to get monkSDF{1}
                clipSpikes = spikeMat;
%                 clipSpikes(spikeMat > repmat(Task.SRT+Task.GoCue,1,size(spikeMat,2))) = nan;
                if sum(isfinite(clipSpikes(:))) < 5
                    monkSDF{ir,1} = nan(1,length(vWind));
                    monkSDFTimes{ir,1} = vWind;
                else
                    [cSDF,cTimes] = klSpkRatev2(clipSpikes,'-q',1,'-k',kType,'-w',gWidth);
                    monkSDF{ir,1} = nanmean(cSDF(Task.TargetLoc==rf(1),ismember(cTimes,vWind)),1);
                    monkSDFTimes{ir,1} = cTimes(ismember(cTimes,vWind));
%                     monkSDF{ir,1} = nanmean(vSDF(Task.TargetLoc==rf(1),ismember(vTimes,vWind)),1);
%                     monkSDFTimes{ir,1} = vTimes(ismember(vTimes,vWind));
                end
            case {'Search','search'}
                % Row 1: Target in RF
                % Row 2: Non-Salient Distractor in RF
                % Row 3: Salient Distractor in RF
                % Row 4: Any Distractor in RF
                monkSDF{ir,2} = nan(2,length(mWind));
                monkSDF{ir,1} = nan(2,length(vWind));
                if ~badM
                    if any(Task.TargetLoc==rf(2) & Task.Correct == 1 & Task.Singleton==0)
                        monkSDF{ir,2}(1,:) = nanmean(mSDF(Task.TargetLoc==rf(2) & Task.Correct == 1 & Task.Singleton==0,ismember(mTimes,mWind)),1);
                    end
                    %monkSDF{ir,2}(2,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc==rf(2) & Task.Correct == 1,ismember(mTimes,mWind)),1);
                    %monkSDF{ir,2}(3,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc~=rf(2) & Task.Correct == 1,ismember(mTimes,mWind)),1);
                    if any(Task.TargetLoc ~= rf(2) & Task.Correct==1 & Task.Singleton==0)
                        monkSDF{ir,2}(2,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.Correct == 1 & Task.Singleton==0,ismember(mTimes,mWind)),1);
                    end
                    monkSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
                else
                    monkSDF{ir,2} = nan(4,length(mWind));
                    monkSDFTimes{ir,2} = mWind;
                end
                if ~badV
                    if any(Task.TargetLoc==rf(1) & Task.Correct==1 & Task.Singleton==0)
                        monkSDF{ir,1}(1,1:sum(ismember(vTimes,vWind))) = nanmean(vSDF(Task.TargetLoc==rf(1) & Task.Correct == 1 & Task.Singleton==0,ismember(vTimes,vWind)),1);
                    end
                    %monkSDF{ir,1}(2,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.DistLoc==rf(1) & Task.Correct == 1,ismember(vTimes,vWind)),1);
                    %monkSDF{ir,1}(3,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.DistLoc~=rf(1) & Task.Correct == 1,ismember(vTimes,vWind)),1);
                    if any(Task.TargetLoc~=rf(1) & Task.Correct==1 & Task.Singleton==0)
                        monkSDF{ir,1}(2,1:sum(ismember(vTimes,vWind))) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.Correct == 1 & Task.Singleton==0,ismember(vTimes,vWind)),1);
                    end
                    monkSDFTimes{ir,1} = vWind;%vTimes(ismember(vTimes,vWind));
                else
                    monkSDF{ir,1} = nan(4,length(vWind));
                    monkSDFTimes{ir,1} = vWind;
                end
                if doReward
                    if ~badR
                        monkSDF{ir,3}(1,:) = nanmean(rSDF(Task.TargetLoc==rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                        monkSDF{ir,3}(2,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.DistLoc==rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                        monkSDF{ir,3}(3,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.DistLoc~=rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                        monkSDF{ir,3}(4,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                        monkSDFTimes{ir,3} = rTimes(ismember(rTimes,rWind));
                    else
                        monkSDF{ir,3} = nan(4,length(rWind));
                        monkSDFTimes{ir,3} = rWind;
                    end
                end
            case {'Capture','capture','Cap','cap'}
                % Row 1: Target in RF, no salient distractor
                % Row 2: Distractor in RF, no salient distractors
                % Row 3: Target in RF, salient distractor
                % Row 4: Non-salient distractor in RF, salient distractor
                % present on trial
                % Row 5: Salient distractor in RF

                if ~badM
                    monkSDF{ir,2}(1,:) = nanmean(mSDF(Task.TargetLoc==rf(2) & Task.Correct==1 & Task.Singleton==0,ismember(mTimes,mWind)),1);
                    monkSDF{ir,2}(2,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.Correct==1 & Task.Singleton==0,ismember(mTimes,mWind)),1);
                    monkSDF{ir,2}(3,:) = nanmean(mSDF(Task.TargetLoc==rf(2) & Task.Correct==1 & Task.Singleton==1,ismember(mTimes,mWind)),1);
                    monkSDF{ir,2}(4,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc ~= rf(2) & Task.Correct==1 & Task.Singleton==1,ismember(mTimes,mWind)),1);
                    monkSDF{ir,2}(5,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc == rf(2) & Task.Correct==1 & Task.Singleton==1,ismember(mTimes,mWind)),1);
                    monkSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
                else
                    monkSDF{ir,2} = nan(5,length(mWind));
                    monkSDFTimes{ir,2} = mWind;
                end
                if ~badV
                    monkSDF{ir,1}(1,:) = nanmean(vSDF(Task.TargetLoc==rf(1) & Task.Correct==1 & Task.Singleton==0,ismember(vTimes,vWind)),1);
                    monkSDF{ir,1}(2,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.Correct==1 & Task.Singleton==0,ismember(vTimes,vWind)),1);
                    monkSDF{ir,1}(3,:) = nanmean(vSDF(Task.TargetLoc==rf(1) & Task.Correct==1 & Task.Singleton==1,ismember(vTimes,vWind)),1);
                    monkSDF{ir,1}(4,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.DistLoc ~= rf(1) & Task.Correct==1 & Task.Singleton==1,ismember(vTimes,vWind)),1);
                    monkSDF{ir,1}(5,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.DistLoc == rf(1) & Task.Correct==1 & Task.Singleton==1,ismember(vTimes,vWind)),1);
                    monkSDFTimes{ir,1} = vTimes(ismember(vTimes,vWind));
                else
                    monkSDF{ir,1} = nan(5,length(vWind));
                    monkSDFTimes{ir,1} = vWind;
                end
                if doReward
                    if ~badR
                        monkSDF{ir,3}(1,:) = nanmean(rSDF(Task.TargetLoc==rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                        monkSDF{ir,3}(2,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.DistLoc==rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                        monkSDF{ir,3}(3,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.DistLoc~=rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                        monkSDF{ir,3}(4,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                        monkSDFTimes{ir,3} = rTimes(ismember(rTimes,rWind));
                    else
                        monkSDF{ir,3} = nan(4,length(rWind));
                        monkSDFTimes{ir,3} = rWind;
                    end
                end
        end
        fprintf(repmat('\b',1,length(printStr)+length(printHead)));
    end
    allSDF = cat(1,allSDF,monkSDF);
    allSDFTimes = cat(1,allSDFTimes,monkSDFTimes);
    allAreas = cat(1,allAreas,monkAreas);
    allMonks = cat(1,allMonks,monkMonks);
    allRows = cat(1,allRows,rowInds);
    allOtherTasks = cat(1,allOtherTasks,monkOtherTasks);
    allSess = cat(1,allSess,monkSess);
    allChans = cat(1,allChans,monkChans);
    if ~exist('allParams','var')
        allParams = monkParams;
    else
        allParams((1:length(monkParams))+length(allParams)) = monkParams;
    end
end

end

function gFit = fitGauss(Task,sdf,times,wind)

uLocs = nunique(Task.TargetLoc);
% Cut locations that don't make sense.
[n,c] = hist(mod(uLocs,45),nunique(mod(uLocs,45)));
uLocs(mod(uLocs,45)~=c(n==max(n))) = [];
resp = nan(size(uLocs));
for il = 1:length(uLocs)
    resp(il) = nanmean(nanmean(sdf(Task.TargetLoc==uLocs(il) & Task.Correct==1,ismember(times,wind)),2),1);
end
% Get maximum value
maxLocInd = find(resp==max(resp),1);
shiftAmnt = uLocs(maxLocInd)-180;
blInds = mod((-1:1)+mod(maxLocInd+length(resp)/2,length(resp)),length(resp)); blInds(blInds==0) = length(resp);
blSub = nanmean(resp(blInds));
fitX = mod(uLocs(~isnan(resp)),360);
fitY = resp(~isnan(resp))-blSub;
if sum(isfinite(fitY)) < 3
    [gFit.amp,gFit.mu,gFit.sig,gFit.bl] = deal(nan);
    gFit.gof = struct('sse',nan,'rsquare',nan,'dfe',nan,'adjrsquare',nan,'rmse',nan);
    return
end
try
    [f, gof]=fit(fitX,fitY,'gauss1');
catch
    [f.a1, f.b1, f.c1] = deal(nan);
    gof = struct('sse',nan,'rsquare',nan,'dfe',nan,'adjrsquare',nan,'rmse',nan);
end
gFit.amp = f.a1;
gFit.mu = f.b1+shiftAmnt;
gFit.sig = f.c1;
gFit.bl = blSub;
gFit.gof = gof;
end