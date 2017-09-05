function rateStruct = klRateStats(spiketimes)

% Set defaults
blWind      = [-300,-100];

% Set filters for nonsense values
rateLim     = 100;
visMax      = 200;
movMin      = -200;
alph        = .1;

%% Get rate statistics for all spikes
[vSDF,vTimes] = klSpkRatev2(spiketimes,'-q',1);
blMean = nanmean(nanmean(vSDF(:,vTimes >= blWind(1) & vTimes < blWind(2)),1),2);
blStd = nanstd(nanmean(vSDF(:,vTimes >= blWind(1) & vTimes < blWind(2)),1),[],2);
            
[mnRate, stdRate] = klGetMeanRate(spiketimes);
if blMean > rateLim,
    mnRate = nan; stdRate = nan; blMean = nan;
end

%             fano              = SPK_get_Fano(spiketimes);
fano              = klGetFano(spiketimes);
cv                = klGetCV(spiketimes);
cv2               = klGetCV(spiketimes,'-type','local');
lv                = klGetLV(spiketimes);
lvr               = klGetLV(spiketimes,'-type','revised');

% Get ISI stats
%             blSpikes          = spiketimes; blSpikes(blSpikes > 0) = nan;
isiMat            = diff(spiketimes,1,2);
meanISI           = nanmean(isiMat(:));
stdISI            = nanstd(isiMat(:));
isiShort          = sum(isiMat(:) < 3)/sum(isfinite(isiMat(:)));

rateStruct.trial.mnRate = mnRate;
rateStruct.trial.stdRate = stdRate;
rateStruct.trial.blMean = blMean;
rateStruct.trial.blStd = blStd;
rateStruct.trial.fano = fano;
rateStruct.trial.cv = cv;
rateStruct.trial.cv2 = cv2;
rateStruct.trial.lv = lv;
rateStruct.trial.lvr = lvr;
rateStruct.trial.mnISI = meanISI;
rateStruct.trial.stdISI = stdISI;
rateStruct.trial.isiShort = isiShort;


%% Now for baseline spikes
blSpks = spiketimes;
blSpks(blSpks < blWind(1) | blSpks > blWind(2)) = nan;

[mnRate, stdRate] = klGetMeanRate(blSpks);
if blMean > rateLim,
    mnRate = nan; stdRate = nan; blMean = nan;
end

%             fano              = SPK_get_Fano(spiketimes);
fano              = klGetFano(blSpks);
cv                = klGetCV(blSpks);
cv2               = klGetCV(blSpks,'-type','local');
lv                = klGetLV(blSpks);
lvr               = klGetLV(blSpks,'-type','revised');

% Get ISI stats
%             blSpikes          = spiketimes; blSpikes(blSpikes > 0) = nan;
isiMat            = diff(blSpks,1,2);
meanISI           = nanmean(isiMat(:));
stdISI            = nanstd(isiMat(:));
isiShort          = sum(isiMat(:) < 3)/sum(isfinite(isiMat(:)));

rateStruct.baseline.mnRate = mnRate;
rateStruct.baseline.stdRate = stdRate;
rateStruct.baseline.blMean = blMean;
rateStruct.baseline.blStd = blStd;
rateStruct.baseline.fano = fano;
rateStruct.baseline.cv = cv;
rateStruct.baseline.cv2 = cv2;
rateStruct.baseline.lv = lv;
rateStruct.baseline.lvr = lvr;
rateStruct.baseline.mnISI = meanISI;
rateStruct.baseline.stdISI = stdISI;
rateStruct.baseline.isiShort = isiShort;
