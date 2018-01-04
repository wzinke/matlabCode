function outCell = klProAntiPhysioStats(sessName,varargin)

baseDir = '/mnt/teba/Users/Kaleb/proAntiProcessed/';
doSave = 1;
watch = 0;
fresh = 0;
subSamp = 1;
consecCrit = 20;
alph = .05;
makePlots = 0;
isiCheck = 2;
doStats = 1;
writeOut = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-f'}
            fresh = varargin{varStrInd(iv)+1};
    end
end

% Set up header
header = {'Session','ChanID','Channel','Unit',...
    'SNR','WvWidth',',TfR','mnRate','mnBaseline','ISI<2ms',...
    'SSTp','SSTa','ESTa','SRT',...
    'vTT','vTL','vCorr','vE','vCong','TTxvTL','vTTxvCorr','vTTxvE','vTTxvCong','vTLxvCorr','vTLxvE','vTLxvCong','vCorrxvE','vCorrxvCong','vExvCong','mTT','mCorr','mE','mCong','mEnd','mTTxmCorr','mTTxmE','mTTxmCong','mTTxmEnd','mCorrxmE','mCorrxmCong','mCorrxmEnd','mExmCong','mExmEnd','mCongxmEnd'};
statsFile = './proAntiBookKeepingv2.xls';

% Load in behavior and eye positions
sessDir = [baseDir,sessName];
load([sessDir, '/Behav.mat']);
load([sessDir, '/Eyes.mat']);

% Find unit .mat files
chanDirs = dir([sessDir,'/Chan*']);

% Set window for visual RF determination
vWind = 50:150;
blWind = -300:-200;
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

% Get unique target locations
proLocs = nunique(Task.TargetLoc(isPA));
antiLocs = mod(proLocs+180,360);
uEccs = nunique(Task.Eccentricity(isPA));
if length(proLocs)==4
    spInds = [6,2,4,8];
end

cutTask = klCutTask(Task,isPA & isGood);
vVars = {'TrialType','TargetLoc','Correct','Eccentricity','Congruent'};
vFactMat = cell(1,length(vVars));
for iv = 1:length(vVars)
    vFactMat{iv} = cutTask.(vVars{iv});
end
mVars = {'TrialType','Correct','Eccentricity','Congruent','EndStimInd'};
mFactMat = cell(1,length(mVars));
for im = 1:length(mVars)
    mFactMat{im} = cutTask.(mVars{im});
end

outCell = {};
for ic = 1:length(chanDirs)
    fprintf('Analyzing Channel %d (of %d)...\n',ic,length(chanDirs));
    clearvars -except ic chanMats chanDirs Task Eyes sessDir sessName... 
        is* proLocs antiLocs congMat congStyle spInds uEccs ...
        vWind vXwind mXwind vVars mVars vFactMat mFactMat blWind...
        doSave watch fresh makePlots...
        subSamp consecCrit...
        doStats alph...
        outCell statsFile header writeOut
        
        
    % Get some identifying info
    chanMats = dir([sessDir,filesep,chanDirs(ic).name,'/chan*.mat']);
    chanCode = chanMats.name;
    channel = chanCode((strfind(chanCode,'chan')+4):(strfind(chanCode,'.mat')-2));
    unit = abc2num(chanCode(strfind(chanCode,'.mat')-1));
    
    % Load in spikes
    load([chanMats.folder,filesep,chanMats.name]);
    spkMat = klPlaceEvents(Task,spkTimes);
    
    % Get some waveform info
    [snr] = klGetSNRv1(waves);
    [width,tfr] = klWvWidthv3(nanmean(waves,1),(-25:25).*(1e6/24414));
    
    % Get the ISI "distribution"
    isiVect = diff(spkTimes);
    percISI = sum(isiVect <= isiCheck)./sum(isfinite(isiVect));
    
    % If not enough spikes, move on
    if nanmean(nansum(isfinite(spkMat),2),1) < 2.5
        continue
    end
    
    % Make SDFs
    [vSDF,vTimes] = klSpkRatev2(spkMat,'-q',1);
    [mSDF,mTimes] = klSpkRatev2(spkMat-repmat(Task.SRT,1,size(spkMat,2)),'-q',1);
    [rSDF,rTimes] = klSpkRatev2(spkMat-repmat(Task.Reward,1,size(spkMat,2)),'-q',1);
    [tSDF,tTimes] = klSpkRatev2(spkMat-repmat(Task.Tone,1,size(spkMat,2)),'-q',1);
    
    % Get receptive field from pro trials
    for il = 1:length(proLocs)
        vResp(il) = nanmean(nanmean(vSDF(isPro & isCorrect & isPA & Task.TargetLoc==proLocs(il),ismember(vTimes,vWind))));
    end
    [maxVals,maxInd] = max(vResp);
    oppInd = mod(maxInd+length(proLocs)/2,length(proLocs)); if oppInd == 0, oppInd = length(proLocs); end
    vRF = proLocs(maxInd);
    
    
    % Cut down vSDF and mSDF to the times of interest
    vSub = vTimes(ismember(vTimes,[vXwind(1):vXwind(2)]));
    mSub = mTimes(ismember(mTimes,[mXwind(1):mXwind(2)]));
    vCut = vSDF(isPA & isGood,ismember(vTimes,[vXwind(1):vXwind(2)]));
    mCut = mSDF(isPA & isGood,ismember(mTimes,[mXwind(1):mXwind(2)]));
    
    % Get mean firing rates
    blMean = nanmean(nanmean(vSDF(isPA & isGood,ismember(vTimes,blWind)),1),2);
    mnRate = nanmean(nanmean([vCut,mCut],1),2);
    
    % Run some stats on relevant factors
    if doStats
        for ii = 1:subSamp:size(vCut,2)
            [pv(:,ii)] = anovan(vCut(:,ii),vFactMat,'model','interaction','varnames',vVars,'display','off');
            [pm(:,ii)] = anovan(mCut(:,ii),mFactMat,'model','interaction','varnames',mVars,'display','off');
        end

        % Find times of significant effects (via ANOVA)
        consecV = klGetConsecutive(pv < alph);
        consecM = klGetConsecutive(pm < alph);
        vMeets = nan(1,size(pv,1));
        mMeets = nan(1,size(pm,1));
        for iv = 1:size(pv,1)
            tmp = find(consecV > (consecCrit/subSamp),1);
            if isempty(tmp)
                vMeets(iv) = nan;
            else
                vMeets(iv) = tmp;
            end
        end
        for im = 1:size(pm,1)
            tmp = find(consecM > (consecCrit/subSamp),1);
            if isempty(tmp)
                mMeets(iv) = nan;
            else
                mMeets(iv) = tmp;
            end
        end
    else
        vMeets = nan(sum(1:size(vFactMat,2)),1);
        mMeets = nan(sum(1:size(mFactMat,2)),1);
    end
    
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
    
    
    
    % Put everything in the output
    outRow = {sessName,chanCode,channel,unit,snr,width,tfr,mnRate,blMean,percISI,sstp,ssta,esta,srt};
    for iv = 1:length(vMeets)
        outRow{length(outRow)+1} = vMeets(iv);
    end
    for im = 1:length(mMeets)
        outRow{length(outRow)+1} = mMeets(im);
    end
    outCell = cat(1,outCell,outRow);
    
    % Send info to make plots if asked
    if makePlots
        klProAntiPhysio(sessName,'-s',doSave);
        if watch
            keyboard
        end
        close all
    end
end

% Write output
if writeOut
    if exist(statsFile,'file') && ~fresh
        [~,~,old] = xlsread(statsFile);
        outCell = cat(1,old,outCell);
    else
        outCell = cat(1,header,outCell);
    end
    xlwrite(statsFile,outCell);    
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
    
    
    
    
    
    
    
    
    
    
    
    
    