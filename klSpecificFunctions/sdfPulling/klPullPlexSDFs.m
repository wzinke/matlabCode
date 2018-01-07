function [allSDF,allSDFTimes,allAreas] = klPullPlexSDFs(task,varargin)


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

% Set constants
minRate = 5;
isiCrit = inf;
snrCrit = .85;
areaCrit = {'FEF','F2'};

% Decode varargin

% Set helpful windows
vWind = -300:500;
mWind = -500:300;

monks = {'Gauss','Helmholtz','Darwin'};
allSDF = {};
allSDFTimes = {};
allAreas = {};
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
    goodRows = cell2mat(all(5:end,snrCol)) > snrCrit & cell2mat(all(5:end,mnRateCol)) > minRate & cell2mat(all(5:end,isiColRate)) < isiCrit & ismember(all(5:end,areaCol),areaCrit);
    rowInds = find(goodRows)+headRow;
    
    monkSDF = cell(length(rowInds),3);
    monkSDFTimes = cell(length(rowInds),3);
    monkAreas = all(rowInds,areaCol);
    % Start row loop
    for ir = 1:length(rowInds)
        myRow = rowInds(ir);
        printHead = sprintf('Analyzing row %d (%d of %d):  ',myRow,ir,length(rowInds));
        fprintf(printHead);

        % Find and load the appropriate .mat files for spikes and behavior
        printStr = 'Loading data...';
        fprintf(printStr);
        
        switch task
            case {'MG','mg'}
                fileName = all{myRow,4};
            case {'Search','Capture','search','capture','Cap','cap'}
                fileTmp = dir(sprintf('%s/%s/DSP/%s/*Cap*.mat',monkFolder,all{myRow,sessCol},all{myRow,chanCol}));
                if isempty(fileTmp)
                    fileTmp = dir(sprintf('%s/%s/DSP/%s/*Search*.mat',monkFolder,all{myRow,sessCol},all{myRow,chanCol}));
                end
                if isempty(fileTmp)
                    continue
                else
                    fileName = fileTmp.name;
                end
        end
        fileStr = sprintf('%s/%s/DSP/%s/%s',monkFolder,all{myRow,sessCol},all{myRow,chanCol},fileName);
        load(fileStr);
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
            case {'Search','Capture','search','capture','Cap','cap'}
                goodTrials = ismember(Task.TaskType,{'Search','Cap'});
        end
        spikeMat = spikeMat(goodTrials,:);
        Task = klCutTask(Task,goodTrials);
        % Get SDFs
        badV = 0; badM = 0;
        if sum(isfinite(spikeMat(:))) < 5
            vSDF = 0; vTimes = 0; badV = 1;
        else
            [vSDF,vTimes] = klSpkRatev2(spikeMat,'-q',1);
        end
        mMat = spikeMat-repmat(Task.SRT + Task.GoCue,1,size(spikeMat,2));
        if sum(isfinite(mMat(:))) < 5
            mSDF = 0; mTimes = 0; badM = 1;
        else
            [mSDF,mTimes] = klSpkRatev2(spikeMat-repmat(Task.SRT + Task.GoCue,1,size(spikeMat,2)),'-q',1);
        end

        % Get RF
        rf = getRFv2(spikeMat,Task,'-s',1);
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
                if ~badV
                    monkSDF{ir,3} = nanmean(vSDF(Task.TargetLoc==rf(1),ismember(vTimes,vWind)),1);
                    monkSDFTimes{ir,3} = vTimes(ismember(vTimes,vWind));
                else
                    monkSDF{ir,3} = nan(1,length(vWind));
                    monkSDFTimes{ir,3} = vWind;
                end

                % Clip post-mov spikes from vis to get monkSDF{1}
                clipSpikes = spikeMat;
                clipSpikes(spikeMat > repmat(Task.SRT+Task.GoCue,1,size(spikeMat,2))) = nan;
                if sum(isfinite(clipSpikes(:))) < 5
                    monkSDF{ir,1} = nan(1,length(vWind));
                    monkSDFTimes{ir,1} = vWind;
                else
                    [cSDF,cTimes] = klSpkRatev2(clipSpikes,'-q',1);
                    monkSDF{ir,1} = nanmean(cSDF(Task.TargetLoc==rf(1),ismember(cTimes,vWind)),1);
                    monkSDFTimes{ir,1} = cTimes(ismember(cTimes,vWind));
                end
            case {'Search','Capture','search','capture','Cap','cap'}
                % Row 1: Target in RF
                % Row 2: Non-Salient Distractor in RF
                % Row 3: Salient Distractor in RF
                % Row 4: Any Distractor in RF
                if ~badM
                    monkSDF{ir,2}(1,:) = nanmean(mSDF(Task.TargetLoc==rf(2) & Task.Correct == 1,ismember(mTimes,mWind)),1);
                    monkSDF{ir,2}(2,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc==rf(2) & Task.Correct == 1,ismember(mTimes,mWind)),1);
                    monkSDF{ir,2}(3,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc~=rf(2) & Task.Correct == 1,ismember(mTimes,mWind)),1);
                    monkSDF{ir,2}(4,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.Correct == 1,ismember(mTimes,mWind)),1);
                    monkSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
                else
                    monkSDF{ir,2} = nan(4,length(mWind));
                    monkSDFTimes{ir,2} = mWind;
                end
                if ~badV
                    monkSDF{ir,1}(1,:) = nanmean(vSDF(Task.TargetLoc==rf(2) & Task.Correct == 1,ismember(vTimes,vWind)),1);
                    monkSDF{ir,1}(2,:) = nanmean(vSDF(Task.TargetLoc~=rf(2) & Task.DistLoc==rf(2) & Task.Correct == 1,ismember(vTimes,vWind)),1);
                    monkSDF{ir,1}(3,:) = nanmean(vSDF(Task.TargetLoc~=rf(2) & Task.DistLoc~=rf(2) & Task.Correct == 1,ismember(vTimes,vWind)),1);
                    monkSDF{ir,1}(4,:) = nanmean(vSDF(Task.TargetLoc~=rf(2) & Task.Correct == 1,ismember(vTimes,vWind)),1);
                    monkSDFTimes{ir,1} = vTimes(ismember(vTimes,vWind));
                else
                    monkSDF{ir,1} = nan(4,length(mWind));
                    monkSDFTimes{ir,1} = mWind;
                end
        end
        fprintf(repmat('\b',1,length(printStr)+length(printHead)));
    end
    allSDF = cat(1,allSDF,monkSDF);
    allSDFTimes = cat(1,allSDFTimes,monkSDFTimes);
    allAreas = cat(1,allAreas,monkAreas);
end

