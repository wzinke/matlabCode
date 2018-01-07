function [outSDF,outSDFTimes,outAreas,outRows,outMonks,outOtherTasks,outSess,outChans,outParams,outVis,outTST,outDST,outPerf,outRF,outRT,outTSTA] = klPullKilosortSDFs(task,varargin)

% Set data location
baseDir = '/media/loweka/ExtraDrive1/dataProcessed/';

% Set constants
minRate = 5;
isiCrit = inf;
snrCrit = 1;
areaCrit = {'FEF','F2'};
doReward = 0;
kType = 'psp';
gWidth = 20;
clip = 1;
waveRate = 24414;

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
        case {'rows'}
            manRows = varargin{varStrInd(iv)+1};
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

% Load in excel file
[num,text,all]=xlsread('kiloSortBookKeeping.xlsx');

% Pick out good columns
snrs = cellfun(@(x) single(x), all(:,3),'UniformOutput',0);
mnRates = cellfun(@(x) single(x), all(:,5),'UniformOutput',0);
goodISI = cellfun(@(x) single(x), all(:,6),'UniformOutput',0);

goodRows = cellfun(@(x) x > snrCrit,snrs) & cellfun(@(x) x > minRate,mnRates) & cellfun(@(x) x < isiCrit, goodISI) & ismember(all(:,7),areaCrit);
rowInds = find(goodRows);
if exist('manRows','var')
    rowInds = manRows;
end
if doReward
    outSDF = cell(length(rowInds),3);   
    outSDFTimes = cell(length(rowInds),3);
else
    outSDF = cell(length(rowInds),2);   
    outSDFTimes = cell(length(rowInds),2);
end

outAreas = all(rowInds,7);
outMonks = cell(length(rowInds),1); [outMonks{1:length(rowInds)}] = deal('Darwin');
outRows = rowInds;
outVis = false(length(rowInds),1);
outTST = nan(length(rowInds),1);
outDST = nan(length(rowInds),1);
outTSTA = nan(length(rowInds),1);
outPerf = nan(length(rowInds),2);
outRF = nan(length(rowInds),2);
outRT = nan(length(rowInds),2);
outSess = cellfun(@(x) x(2:(end-1)),all(rowInds,1),'UniformOutput',0);
outChans = cellfun(@(x) x(2:(end-5)),all(rowInds,2),'UniformOutput',0);
outOtherTasks = cell(length(rowInds),1);
[outParams(1:length(rowInds))] = deal(struct('theta',[],'sigma',[],'alpha',[],'offset',[],'blRate',[],'Fano',[],'CV',[],'CV2',[],'LV',[],'LVR',[],'vGOF',[],'mGOF',[],'spkWidth',[]));
    
% Loop through rows
for ir = 1:length(rowInds)
    printHead = sprintf('Analyzing row %d of %d:  ',ir,length(rowInds));
    fprintf(printHead);
    
    clear waves mnWave
    
    % Find and load the appropriate .mat files for spikes and behavior
    printStr = 'Loading data...';
    fprintf(printStr);
    load([[tebaMount,'Users/Kaleb/dataProcessed/',all{rowInds(ir),1}(2:(strfind(all{rowInds(ir),1},'_probe')-1)),'/Behav.mat']);
    load(sprintf('%s/%s/%s',baseDir,all{rowInds(ir),1}(2:(end-1)),all{rowInds(ir),2}(2:(end-1))),'spkTimes','waves');
    fprintf(repmat('\b',1,length(printStr)));
    
    waveTimes = (1:size(waves,2)).*(1e6/waveRate);
    mnWave = nanmean(waves); clear waves;
    splT = spline(1:length(waveTimes),waveTimes,1:.1:length(waveTimes));
    splWave = spline(1:length(mnWave),mnWave,1:.1:length(mnWave));
    
    printStr = 'Calculating SDFs...';
    fprintf(printStr);
    % Calculate SDF for MG and in all MG locations
    spikeMat = klPlaceEvents(Task,spkTimes);
    
%     if length(unique(Task.TaskType)) > 1
%         keyboard
%     end
%     
    
    % Chop down Task a little bit...
    mgTrials = ismember(Task.TaskType,'MG');
    switch task
        case {'MG','mg'}
            goodTrials = ismember(Task.TaskType,'MG');
            otherCut = 'MG';
        case {'Search','search'}
            goodTrials = ismember(Task.TaskType,{'Search','Cap'});
            mgTrials = mgTrials & ismember(Task.TargetLoc,0:90:359);
            otherCut = {'Search','Cap','Capture'};
        case {'Capture','capture','Cap','cap'}
            goodTrials = ismember(Task.TaskType,{'Cap'});
            otherCut = {'Cap','Capture'};
            mgTrials = mgTrials & ismember(Task.TargetLoc,0:90:359);
    end
    allTasks = unique(Task.TaskType);
    outOtherTasks{ir} = allTasks(~ismember(allTasks,otherCut));
    hasGoodTrials(ir) = any(goodTrials);
    if ~any(goodTrials)
        fprintf(repmat('\b',1,length(printStr)+length(printHead)));
        continue
    end
    % Get RF
%     if sum(mgTrials) > 50
%         rf = getRFv2(spikeMat(mgTrials,:),klCutTask(Task,mgTrials));
%         Task = klCutTask(Task,goodTrials);
%         spikeMat = spikeMat(goodTrials,:);
%     else
        Task = klCutTask(Task,goodTrials);
        spikeMat = spikeMat(goodTrials,:);
        rf = getRFv2(spikeMat,Task,'-s',1);
%     end
    
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
        clipMat = spikeMat;
        clipMat(clipMat >= repmat(Task.SRT+Task.GoCue,1,size(clipMat,2))) = nan;
        if sum(isfinite(clipMat(:))) < 5
            vSDF = 0; vTimes = 0; badV = 1;
        else
            [vSDF,vTimes] = klSpkRatev2(clipMat,'-q',1,'-k',kType,'-w',gWidth);
        end
    else
        if sum(isfinite(spikeMat(:))) < 5
            vSDF = 0; vTimes = 0; badV = 1;
        else
            [vSDF,vTimes] = klSpkRatev2(spikeMat,'-q',1,'-k',kType,'-w',gWidth);
        end
    end
    fprintf(repmat('\b',1,length(printStr)));
    printStr = 'Putting into output cell...';
    fprintf(printStr);
    
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
    outParams(ir).theta = [vFit.mu,mFit.mu];
    outParams(ir).sigma = [vFit.sig,mFit.sig];
    outParams(ir).alpha = [vFit.amp,mFit.amp];
    outParams(ir).offset = [vFit.bl,mFit.bl];
    outParams(ir).vGOF = vFit.gof;
    outParams(ir).mGOF = mFit.gof;
    outParams(ir).blRate = nanmean(nanmean(vSDF(:,ismember(vTimes,blWind)),2),1);
    outParams(ir).CV = klGetCV(spikeMat);
    outParams(ir).CV2 = klGetCV(spikeMat,'-type','local');
    outParams(ir).LV = klGetLV(spikeMat);
    outParams(ir).LVR = klGetLV(spikeMat,'-type','revised');
    outParams(ir).Fano = klGetFano(spikeMat);
    outParams(ir).spkWidth = abs(splT(splWave==max(splWave))-splT(splWave==min(splWave)));%klWvWidthv3();
    
    %% Get relevant task features
    switch task
        case {'MG','mg'}
            % Pull mov and raw visual SDF into outSDF and outSDFTimes. Just
            % in RF
            if ~badM
                outSDF{ir,2} = nanmean(mSDF(Task.TargetLoc==rf(2),ismember(mTimes,mWind)),1);
                outSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
            else
                outSDF{ir,2} = nan(1,length(mWind));
                outSDFTimes{ir,2} = mWind;
            end
%             if ~badV
%                 outSDF{ir,3} = nanmean(vSDF(Task.TargetLoc==rf(1),ismember(vTimes,vWind)),1);
%                 outSDFTimes{ir,3} = vTimes(ismember(vTimes,vWind));
%             else
%                 outSDF{ir,3} = nan(1,length(vWind));
%                 outSDFTimes{ir,3} = vWind;
%             end
            if doReward
                if ~badR
                    outSDF{ir,3} = nanmean(rSDF(ismember(Task.TargetLoc,rf),ismember(rTimes,rWind)),1);
                    outSDFTimes{ir,3} = rTimes(ismember(rTimes,rWind));
                else
                    outSDF{ir,3} = nan(1,length(rWind));
                    outSDFTimes{ir,3} = rWind;
                end
            end

            % Clip post-mov spikes from vis to get outSDF{1}
%             clipSpikes = spikeMat;
%             clipSpikes(spikeMat > repmat(Task.SRT+Task.GoCue,1,size(spikeMat,2))) = nan;
            if sum(isfinite(spikeMat(:))) < 5
                outSDF{ir,1} = nan(1,length(vWind));
                outSDFTimes{ir,1} = vWind;
            else
%                 [cSDF,cTimes] = klSpkRatev2(clipSpikes,'-q',1,'-k',kType,'-w',gWidth);
                outSDF{ir,1} = nanmean(vSDF(Task.TargetLoc==rf(1),ismember(vTimes,vWind)),1);
                outSDFTimes{ir,1} = vTimes(ismember(vTimes,vWind));
            end
        case {'Search','search'}
            % Row 1: Target in RF
            % Row 2: Non-Salient Distractor in RF
            % Row 3: Salient Distractor in RF
            % Row 4: Any Distractor in RF
            outSDF{ir,2} = nan(2,length(mWind));
            outSDF{ir,1} = nan(2,length(vWind));
            if ~badM
                outSDF{ir,2}(1,:) = nanmean(mSDF(Task.TargetLoc==rf(2) & Task.Correct == 1 & Task.Singleton==0,ismember(mTimes,mWind)),1);
                %outSDF{ir,2}(2,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc==rf(2) & Task.Correct == 1,ismember(mTimes,mWind)),1);
                %outSDF{ir,2}(3,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc~=rf(2) & Task.Correct == 1,ismember(mTimes,mWind)),1);
                outSDF{ir,2}(2,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.Correct == 1 & Task.Singleton==0,ismember(mTimes,mWind)),1);
                outSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
            else
                outSDF{ir,2} = nan(4,length(mWind));
                outSDFTimes{ir,2} = mWind;
            end
            if ~badV
                outSDF{ir,1}(1,1:sum(ismember(vTimes,vWind))) = nanmean(vSDF(Task.TargetLoc==rf(1) & Task.Correct == 1 & Task.Singleton==0,ismember(vTimes,vWind)),1);
                %outSDF{ir,1}(2,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.DistLoc==rf(1) & Task.Correct == 1,ismember(vTimes,vWind)),1);
                %outSDF{ir,1}(3,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.DistLoc~=rf(1) & Task.Correct == 1,ismember(vTimes,vWind)),1);
                outSDF{ir,1}(2,1:sum(ismember(vTimes,vWind))) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.Correct == 1 & Task.Singleton==0,ismember(vTimes,vWind)),1);
                outSDFTimes{ir,1} = vWind;%vTimes(ismember(vTimes,vWind));
            else
                outSDF{ir,1} = nan(4,length(mWind));
                outSDFTimes{ir,1} = mWind;
            end
            if doReward
                if ~badR
                    outSDF{ir,3}(1,:) = nanmean(rSDF(Task.TargetLoc==rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                    outSDF{ir,3}(2,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.DistLoc==rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                    outSDF{ir,3}(3,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.DistLoc~=rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                    outSDF{ir,3}(4,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                    outSDFTimes{ir,3} = rTimes(ismember(rTimes,rWind));
                else
                    outSDF{ir,3} = nan(4,length(rWind));
                    outSDFTimes{ir,3} = rWind;
                end
            end
        case {'Capture','capture','Cap','cap'}
            % Get RF
            rf = getRFv2(spikeMat,Task,'-s',1);
            outRF(ir,:) = rf;

            % Check whether the cell is "visual"
%             blSDF = nanmean(vSDF(:,vTimes >= -50 & vTimes < 0),2);
%             tOnSDF = nanmean(vSDF(:,vTimes >= 50 & vTimes <= 150),2);
%             if sum(isfinite(blSDF))<=5 || sum(isfinite(tOnSDF))<=5
%                 visP = nan;
%             else
%                 visP = ranksum(blSDF,tOnSDF,'tail','left');
%             end
%             outVis(ir) = visP < .01 & nanmean(tOnSDF) > nanmean(blSDF);
            mnSDF = nanmean(vSDF,1);
            blMn = nanmean(mnSDF(vTimes >= -200 & vTimes <= -100),2);
            blStd = nanstd(mnSDF(vTimes >= -200 & vTimes <= -100),[],2);
            zSDF = (mnSDF-blMn)./blStd;
            outVis(ir) = any(klGetConsecutive(zSDF(vTimes >= 0 & vTimes <= 100) >= 5) >= 10);
            
            % Loop through time (start at -50) to find TST and DST
            testSDF = vSDF(:,vTimes >= -50 & vTimes <= nanmean(Task.SRT));
            testTimes = vTimes(:,vTimes >= -50 & vTimes <= nanmean(Task.SRT));
            tstP = nan(1,length(testTimes));
            dstP = nan(1,length(testTimes));
            tstA = nan(1,length(testTimes));
            for ii = 1:length(testTimes)
                xt=testSDF(Task.TargetLoc==rf(1) & Task.Correct==1 & Task.Singleton==2222,ii);
                yt=testSDF(Task.TargetLoc~=rf(1) & Task.Correct==1 & Task.Singleton==2222 & Task.DistLoc~=rf(1),ii);
                xa=testSDF(Task.TargetLoc==rf(1) & Task.Correct==1 & Task.Singleton==0,ii);
                ya=testSDF(Task.TargetLoc~=rf(1) & Task.Correct==1 & Task.Singleton==0,ii);
                xd=testSDF(Task.TargetLoc~=rf(1) & Task.Singleton==2222 & Task.DistLoc==rf(1) & Task.Correct==1,ii);
                yd=testSDF(Task.TargetLoc~=rf(1) & Task.Singleton==2222 & Task.DistLoc~=rf(1) & Task.Correct==1,ii);
                if sum(~isnan(xt)) >= 5 && sum(~isnan(yt)) >= 5
                    tstP(ii) = ranksum(xt,yt,'tail','right');
                end
                if sum(~isnan(xd)) >= 5 && sum(~isnan(yd)) >= 5
                    dstP(ii) = ranksum(xd,yd,'tail','left');
                end
                if sum(~isnan(xa)) >= 5 && sum(~isnan(ya)) >= 5
                    tstA(ii) = ranksum(xa,ya,'tail','right');
                end
            end
            tmpT = find(klGetConsecutive(tstP < .01) >= 10,1);
            tmpD = find(klGetConsecutive(dstP < .01) >= 10,1);
            tmpA = find(klGetConsecutive(tstA < .01) >= 10,1);
            if outVis(ir) && ~isempty(tmpT)
                outTST(ir) = testTimes(tmpT);
            else
                outTST(ir) = nan;
            end
            if outVis(ir) && ~isempty(tmpD)
                outDST(ir) = testTimes(tmpD);
            else
                outDST(ir) = nan;
            end
            if outVis(ir) && ~isempty(tmpA)
                outTSTA(ir) = testTimes(tmpA);
            else
                outTSTA(ir) = nan;
            end
            % Row 1: Target in RF, no salient distractor
            % Row 2: Distractor in RF, no salient distractors
            % Row 3: Target in RF, salient distractor
            % Row 4: Non-salient distractor in RF, salient distractor
            % present on trial
            % Row 5: Salient distractor in RF
            
            if ~badM
                outSDF{ir,2}(1,:) = nanmean(mSDF(Task.TargetLoc==rf(2) & Task.Correct==1 & Task.Singleton==0,ismember(mTimes,mWind)),1);
                outSDF{ir,2}(2,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.Correct==1 & Task.Singleton==0,ismember(mTimes,mWind)),1);
                outSDF{ir,2}(3,:) = nanmean(mSDF(Task.TargetLoc==rf(2) & Task.Correct==1 & Task.Singleton==2222,ismember(mTimes,mWind)),1);
                outSDF{ir,2}(4,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc ~= rf(2) & Task.Correct==1 & Task.Singleton==2222,ismember(mTimes,mWind)),1);
                outSDF{ir,2}(5,:) = nanmean(mSDF(Task.TargetLoc~=rf(2) & Task.DistLoc == rf(2) & Task.Correct==1 & Task.Singleton==2222,ismember(mTimes,mWind)),1);
                outSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
            else
                outSDF{ir,2} = nan(5,length(mWind));
                outSDFTimes{ir,2} = mWind;
            end
            if ~badV
                outSDF{ir,1}(1,:) = nanmean(vSDF(Task.TargetLoc==rf(1) & Task.Correct==1 & Task.Singleton==0,ismember(vTimes,vWind)),1);
                outSDF{ir,1}(2,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.Correct==1 & Task.Singleton==0,ismember(vTimes,vWind)),1);
                outSDF{ir,1}(3,:) = nanmean(vSDF(Task.TargetLoc==rf(1) & Task.Correct==1 & Task.Singleton==2222,ismember(vTimes,vWind)),1);
                outSDF{ir,1}(4,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.DistLoc ~= rf(1) & Task.Correct==1 & Task.Singleton==2222,ismember(vTimes,vWind)),1);
                outSDF{ir,1}(5,:) = nanmean(vSDF(Task.TargetLoc~=rf(1) & Task.DistLoc == rf(1) & Task.Correct==1 & Task.Singleton==2222,ismember(vTimes,vWind)),1);
                outSDFTimes{ir,1} = vTimes(ismember(vTimes,vWind));
            else
                outSDF{ir,1} = nan(5,length(vWind));
                outSDFTimes{ir,1} = vWind;
            end
            if doReward
                if ~badR
                    outSDF{ir,3}(1,:) = nanmean(rSDF(Task.TargetLoc==rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                    outSDF{ir,3}(2,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.DistLoc==rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                    outSDF{ir,3}(3,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.DistLoc~=rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                    outSDF{ir,3}(4,:) = nanmean(rSDF(Task.TargetLoc~=rf(1) & Task.Correct == 1,ismember(rTimes,rWind)),1);
                    outSDFTimes{ir,3} = rTimes(ismember(rTimes,rWind));
                else
                    outSDF{ir,3} = nan(4,length(rWind));
                    outSDFTimes{ir,3} = rWind;
                end
            end
            outPerf(ir,1) = nansum(Task.Correct==1 & Task.Singleton==2222)./sum(Task.Abort~=0 & Task.Singleton==2222);
            outPerf(ir,2) = nansum(Task.Correct==1 & Task.Singleton==0)./sum(Task.Abort~=0 & Task.Singleton==0);
            outRT(ir,1) = nanmedian(Task.SRT(Task.Correct==1 & Task.Singleton==2222));
            outRT(ir,2) = nanmedian(Task.SRT(Task.Correct==1 & Task.Singleton==0));
    end
            
    fprintf(repmat('\b',1,length(printStr)+length(printHead)));
end

end

function gFit = fitGauss(Task,sdf,times,wind)

uLocs = nunique(Task.TargetLoc);
% Cut locations that don't make sense.
[n,c] = hist(mod(uLocs,45),nunique(mod(uLocs,45)));
uLocs(mod(uLocs,45)~=c(n==max(n))) = [];
resp = nan(size(uLocs));
for il = 1:length(uLocs)
    resp(il) = nanmean(nanmean(sdf(Task.TargetLoc==uLocs(il),ismember(times,wind)),2),1);
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
    