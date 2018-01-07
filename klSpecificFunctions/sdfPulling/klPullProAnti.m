function [outSDF,outTimes,outTypes,outRTs] = klPullProAnti()

xlFile = 'proAntiBookKeeping2.xlsx';
headRow = 1;

tebaBase = [tebaMount,'Users/Kaleb/proAntiProcessed'];

% set criteria
snrCrit = 1;
isiCrit = .1;
minRate = 5;

maxSSTp = 300;
maxSSTa = 300;
maxESTa = 300;
maxSRT = 300;

vWind = -200:300;
mWind = -300:200;
vCheck = 50:150;
mCheck = -50:0;

% Read in excel file
[num,txt,all] = xlsread(xlFile);

% Define columns
mnRateCol = find(cellfun(@(x) strcmpi(x,'mnRate'),all(headRow,:)),1);
isiCol = find(cellfun(@(x) strcmpi(x,'ISI<2ms'),all(headRow,:)),1);
snrCol = find(cellfun(@(x) strcmpi(x,'SNR'),all(headRow,:)),1);
% areaCol = find(cellfun(@(x) strcmpi(x,'Area'),all(headRow,:)),1);
sessCol = find(cellfun(@(x) strcmpi(x,'Session'),all(headRow,:)),1);
chanCol = find(cellfun(@(x) strcmpi(x,'Channel'),all(headRow,:)),1);
chanMatCol = find(cellfun(@(x) strcmpi(x,'ChanID'),all(headRow,:)),1);
unitCol = find(cellfun(@(x) strcmpi(x,'Unit'),all(headRow,:)),1);
sstpCol = find(cellfun(@(x) strcmpi(x,'SSTp'),all(headRow,:)),1);
sstaCol = find(cellfun(@(x) strcmpi(x,'SSTa'),all(headRow,:)),1);
estaCol = find(cellfun(@(x) strcmpi(x,'ESTa'),all(headRow,:)),1);
srtCol = find(cellfun(@(x) strcmpi(x,'SRT'),all(headRow,:)),1);



% Use criteria to pick out good rows
mnRates = cellfun(@(x) single(x),all((headRow+1):end,mnRateCol),'UniformOutput',0); [mnRates{cellfun(@isempty,mnRates)}] = deal(nan);
isis = cellfun(@(x) single(x),all((headRow+1):end,isiCol),'UniformOutput',0); [isis{cellfun(@isempty,isis)}] = deal(nan);
snrs = cellfun(@(x) single(x),all((headRow+1):end,snrCol),'UniformOutput',0); [snrs{cellfun(@isempty,snrs)}] = deal(nan);
goodRows = cellfun(@(x) x > minRate,mnRates) &...
    cellfun(@(x) x < isiCrit,isis) &...
    cellfun(@(x) x >= snrCrit,snrs);
%     ismember(all((headRow+1):end,areaCol),areaCrit) &...
rowCrit = find(goodRows)+headRow;

% Check cell types
hasSSTp = cellfun(@(x) x <= maxSSTp,all(rowCrit,sstpCol));
hasSSTa = cellfun(@(x) x <= maxSSTa,all(rowCrit,sstaCol));
hasESTa = cellfun(@(x) x <= maxESTa,all(rowCrit,estaCol));
hasSRT = cellfun(@(x) x <= maxSRT,all(rowCrit,srtCol));

isT1 = hasSSTp & hasSSTa;
isT1Flip = hasSSTp & hasSSTa & hasESTa;
isT2 = hasSSTp & ~hasSSTa & hasESTa;
% isPA = hasSRT;

outTypes = [isT1,isT1Flip,isT2,hasSRT];

% Pull in SDFs
outSDF = {};
for ir = 1:length(rowCrit)
    myRow = rowCrit(ir);
    if exist('printStr','var')
        fprintf(repmat('\b',1,length(printStr)));
    end
    printStr = sprintf('Pulling row %d (%d of %d)\n',myRow,ir,length(rowCrit));
    fprintf(printStr);
    
    % Load in behavior
    load([tebaBase,filesep,all{myRow,sessCol},filesep,'Behav.mat']);
    
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

    % Load in spikes
    load([tebaBase,filesep,all{myRow,sessCol},filesep,'Chan',num2str(all{myRow,chanCol}),num2abc(all{myRow,unitCol}),filesep,all{myRow,chanMatCol}]);
    
    % Place in spike matrix
    spkMat = klPlaceEvents(Task,spkTimes);
    movMat = spkMat-repmat(Task.SRT,1,size(spkMat,2));
    
    badV = 0; badM = 0;
    if sum(~isnan(spkMat(:))) < 5
        badV = 1;
    end
    if sum(~isnan(movMat(:))) < 5
        badM = 1;
    end
    
    % Get SDFs
    [vSDF,vTimes] = klSpkRatev2(spkMat,'-q',1);
    [mSDF,mTimes] = klSpkRatev2(movMat,'-q',1);
    
    
    % Calculate RFs
    uLocs = nunique(Task.TargetLoc);
    vResp = nan(1,length(uLocs));
    mResp = nan(1,length(uLocs));
    for il = 1:length(uLocs)
        vResp(il) = nanmean(nanmean(vSDF(isPro & isCorrect & Task.TargetLoc==uLocs(il),ismember(vTimes,vCheck)),1));
        mResp(il) = nanmean(nanmean(mSDF(isPro & isCorrect & Task.TargetLoc==uLocs(il),ismember(mTimes,mCheck)),1));
    end
    if ~any(~isnan(vResp)), vrf = nan; else vrf=uLocs(vResp==max(vResp)); end; if isnan(vrf) || isempty(vrf), badV = 1; end
    if ~any(~isnan(mResp)), mrf = nan; else mrf=uLocs(mResp==max(mResp)); end; if isnan(mrf) || isempty(mrf), badM = 1; end

    % Make output matrix
    if ~badV
        outSDF{ir,1}(1,:) = nanmean(vSDF(isPro & isCorrect & Task.TargetLoc==vrf,ismember(vTimes,vWind)),1);
        outSDF{ir,1}(2,:) = nanmean(vSDF(isPro & isCorrect & Task.TargetLoc==mod(vrf+180,360),ismember(vTimes,vWind)),1);
        outSDF{ir,1}(3,:) = nanmean(vSDF(isAnti & isCorrect & Task.TargetLoc==vrf,ismember(vTimes,vWind)),1);
        outSDF{ir,1}(4,:) = nanmean(vSDF(isAnti & isCorrect & Task.TargetLoc==mod(vrf+180,360),ismember(vTimes,vWind)),1);
    else
        outSDF{ir,1} = nan(4,length(vWind));
    end
    outTimes{ir,1} = vTimes(ismember(vTimes,vWind));
    
    if ~badM
        outSDF{ir,2}(1,:) = nanmean(mSDF(isPro & isCorrect & Task.TargetLoc==mrf,ismember(mTimes,mWind)),1);
        outSDF{ir,2}(2,:) = nanmean(mSDF(isPro & isCorrect & Task.TargetLoc==mod(mrf+180,360),ismember(mTimes,mWind)),1);
        outSDF{ir,2}(3,:) = nanmean(mSDF(isAnti & isCorrect & Task.TargetLoc==mrf,ismember(mTimes,mWind)),1);
        outSDF{ir,2}(4,:) = nanmean(mSDF(isAnti & isCorrect & Task.TargetLoc==mod(mrf+180,360),ismember(mTimes,mWind)),1);
    else
        outSDF{ir,2} = nan(4,length(mWind));
    end
    outTimes{ir,2} = mTimes(ismember(mTimes,mWind));
    
    outRTs{ir,1,1} = Task.SRT(isPro & isCorrect & Task.TargetLoc==vrf);
    outRTs{ir,2,1} = Task.SRT(isPro & isCorrect & Task.TargetLoc==mrf);
    outRTs{ir,3,1} = Task.SRT(isAnti & isCorrect & Task.TargetLoc==vrf);
    outRTs{ir,4,1} = Task.SRT(isAnti & isCorrect & Task.TargetLoc==mrf);
    outRTs{ir,5,1} = Task.SRT(isPro & isCorrect & Task.TargetLoc==mod(vrf+180,360));
    outRTs{ir,6,1} = Task.SRT(isPro & isCorrect & Task.TargetLoc==mod(mrf+180,360));
    outRTs{ir,7,1} = Task.SRT(isAnti & isCorrect & Task.TargetLoc==mod(vrf+180,360));
    outRTs{ir,8,1} = Task.SRT(isAnti & isCorrect & Task.TargetLoc==mod(mrf+180,360));
    
    outRTs{ir,1,2} = Task.SRT(isPro & ~isCorrect & Task.TargetLoc==vrf);
    outRTs{ir,2,2} = Task.SRT(isPro & ~isCorrect & Task.TargetLoc==mrf);
    outRTs{ir,3,2} = Task.SRT(isAnti & ~isCorrect & Task.TargetLoc==vrf);
    outRTs{ir,4,2} = Task.SRT(isAnti & ~isCorrect & Task.TargetLoc==mrf);
    outRTs{ir,5,2} = Task.SRT(isPro & ~isCorrect & Task.TargetLoc==mod(vrf+180,360));
    outRTs{ir,6,2} = Task.SRT(isPro & ~isCorrect & Task.TargetLoc==mod(mrf+180,360));
    outRTs{ir,7,2} = Task.SRT(isAnti & ~isCorrect & Task.TargetLoc==mod(vrf+180,360));
    outRTs{ir,8,2} = Task.SRT(isAnti & ~isCorrect & Task.TargetLoc==mod(mrf+180,360));
    
    
    
%     if ir == 84
%         keyboard
%     end
end
fprintf('\n');
    
    
    
    
    
    
    
    
    
    