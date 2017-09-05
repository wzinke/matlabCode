function [isAP, isNoise, grpMeans, alTimes, outGrp, alWaves, outIDs] = klSortSpikesv4(waveMat,varargin)
warning off

% Set defaults
upSamp = 1;
align = 1;
times         = ((1:32).*25)-(9*25);
tailPerc = 100;
maxWaves = 5000;
type = 'corr';
visualize = 0;
alignVis = 0;
maxK = 6;
pcaType = 'kmeans';
kReps = 20;
forceK = 0;
forceKVal = 2;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'upsamp','-u'},
            upSamp = varargin{varStrInd(iv)+1};
        case {'-a','align'},
            align = varargin{varStrInd(iv)+1};
        case {'-t','times'},
            times = varargin{varStrInd(iv)+1};
        case {'type'}
            type = varargin{varStrInd(iv)+1};
        case {'-v','vis'}
            visualize = varargin{varStrInd(iv)+1};
        case {'maxK'}
            maxK = varargin{varStrInd(iv)+1};
        case {'-k','k'},
            forceKVal = varargin{varStrInd(iv)+1};
        case {'forcek'},
            forceK = varargin{varStrInd(iv)+1};
    end
end

% Upsample, if desired
if upSamp,
    warning off
    if size(waveMat,1) <= maxWaves,
%         diffWaves = [nan(size(waveMat,1),5),diff(waveMat(:,5:20),[],2),nan(size(waveMat,1),length(21:32))];
        diffWaves = [nan(size(waveMat,1),1),diff(waveMat(:,1:30),[],2),nan(size(waveMat,1),length(31:32))];
        cutWaves=waveMat;
        cutWaves(abs(diffWaves) < 2 & ((waveMat > max(waveMat(:)) - (range(waveMat(:))/10)) | waveMat < min(waveMat(:)) + (range(waveMat(:))/10))) = nan;
        for ii = 1:size(waveMat,1),
            smoothWaves(ii,:) = spline(1:32,cutWaves(ii,:),1:.1:32);
        end
        smoothTimes = spline(1:32,times,1:.1:32);
        alignWind = 10;
    else
        if size(waveMat,1) <= (maxWaves),
            [isAP,isNoise, apMean, apTimes] = klSortSpikesv2(waveMat,'-u',0,'-a',1,'-t',times,'type',type);
            spkWaves = waveMat(isAP,:); 
            nzTemp   = waveMat(~isAP,:);
        else
            [sortThresh,sortInd] = sort(abs(nanmean(waveMat(:,times > -50 & times < 50),2)));
            nzTemp = waveMat(sortInd(1:floor(size(waveMat,1)/2)),:);
            spkWaves = waveMat(sortInd((floor(size(waveMat,1)/2)+1):end),:);
        end
        maxSpk = ceil(maxWaves*(size(spkWaves,1)./(size(spkWaves,1)+size(nzTemp,1))));
        maxNz = maxWaves-maxSpk;
        
        spkRand = randperm(size(spkWaves,1));
        spkInds = spkRand(1:min([maxSpk,size(spkWaves,1)]));
        spkClust = spkWaves(spkInds,:);
        nzRand = randperm(size(nzTemp,1));
        nzInds = nzRand(1:min([maxNz,size(nzTemp,1)]));
        nzClust  = nzTemp(nzInds,:);

        subWaves = [spkClust;nzClust];
%         diffWaves = [nan(size(subWaves,1),5),diff(subWaves(:,5:20),[],2),nan(size(subWaves,1),length(21:32))];
        diffWaves = [nan(size(subWaves,1),1),diff(subWaves(:,1:30),[],2),nan(size(subWaves,1),length(31:32))];
        cutWaves=subWaves;
%         cutWaves(abs(diffWaves) < 1 & abs(subWaves) > nanstd(nanmean(subWaves,1))*2.5) = nan;
        cutWaves(abs(diffWaves) < 2 & ((subWaves > max(subWaves(:)) - (range(subWaves(:))/10)) | subWaves < min(subWaves(:)) + (range(subWaves(:))/10))) = nan;
        smoothWaves = nan(size(subWaves,1),length(1:.1:32));
        for ii = 1:size(subWaves,1),
            smoothWaves(ii,:) = spline(1:32,cutWaves(ii,:),1:.1:32);
        end
        smoothTimes = spline(1:32,times,1:.1:32);
        alignWind = 10;
        
        spkClust = smoothWaves(1:size(spkClust,1),:);
        nzClust = smoothWaves((size(spkClust,1)+1):end,:);
%         clear spkWaves nzTemp maxSpk maxNz
%         smoothWaves = [spkClust;nzClust];
    end
    warning on
else
    smoothWaves = waveMat;
    smoothTimes = times;
    alignWind = 2;
end

% if upSamp,
%     diffWaves = [nan(size(waveMat,1),5),diff(waveMat(:,5:20),[],2),nan(size(waveMat,1),length(21:32))];
%     cutWaves=waveMat;
%     cutWaves(abs(diffWaves) < 2 & abs(waveMat) > nanstd(nanmean(waveMat,1))*2.5) = nan;
%     for ii = 1:size(waveMat,1),
%         smoothWaves(ii,:) = spline(1:32,cutWaves(ii,:),1:.1:32);
%     end
%     smoothTimes = spline(1:32,times,1:.1:32);
%     alignWind = 10;
% else
%     smoothWaves = waveMat;
%     smoothTimes = times;
%     alignWind = 4;
% end

% Align these waves on their troughs
if align,
    [alWaves, alTimes] = klTroughAlignv4(smoothWaves,smoothTimes,0,'-w',alignWind,'-v',alignVis);
else
    alWaves = smoothWaves;
    alTimes = smoothTimes;
end

% Get and sort the zero-points
[sortTrough, sortInd] = sort(abs(alWaves(:,alTimes == 0)));

% Get top and bottom n% (or n waves)
if tailPerc <= 1,
    nTake = ceil(tailPerc.*size(alWaves,1));
else
    nTake = round(tailPerc);
end

noiseStart = alWaves(sortInd(1:nTake),:);
signalStart = alWaves(sortInd((end-nTake+1):end),:);
noiseSeed = nanmean(alWaves(sortInd(1:nTake),:),1);
signalSeed = nanmean(alWaves(sortInd((end-nTake+1):end),:),1);
for ii = 1:size(alWaves,2),
    goodCols(ii) = sum(~isnan(alWaves(:,ii))) == size(alWaves,1);
end

switch type
    case {'corr','correlation'},
        
        % Go through each wave and get correlations
        % Column 1 is correlation to signal, column 2 is correlation to noise
        wvCorrs = corr(alWaves(:,goodCols)',[signalSeed(goodCols);noiseSeed(goodCols)]');

        % for iw = 1:size(alWaves,1),
        %     goodInds = ~isnan(noiseSeed) & ~isnan(signalSeed) & ~isnan(alWaves(iw,:));
        %     wvCorrs(iw,1) = corr(alWaves(iw,goodInds)',signalSeed(goodInds)');
        %     wvCorrs(iw,2) = corr(alWaves(iw,goodInds)',noiseSeed(goodInds)');
        % end

        isAP = wvCorrs(:,1) > wvCorrs(:,2);
        isNoise = wvCorrs(:,2) > wvCorrs(:,1);
    case {'pca'}
        for ii = 1:size(smoothWaves,2),
            goodSmooth(ii) = sum(~isnan(smoothWaves(:,ii))) == size(smoothWaves,1);
        end
        % Get PCAs
        [coeffs,scores] = pca(smoothWaves(:,goodSmooth),'NumComponents',3);
        ids = klAgglomv7(scores,'-k',2,'-t','euc');
        [kHat,kGap,kErr,kStruct] = klGetGapv8(scores,ids,'reclust',pcaType,'-t','euc','randreps',10,'-r','rand','randval','unif');
        
        if corr(nanmean(alWaves(ids == 1,goodCols),1)',signalSeed(goodCols)') > corr(nanmean(alWaves(ids == 1,goodCols),1)',noiseSeed(goodCols)'),
            isAP = ids == 1;
            isNoise = ids == 2;
        else
            isAP = ids == 2;
            isNoise = ids == 1;
        end
    case {'kmeans'}
        % Check for all-nan rows
        for ii = 1:size(alWaves,1),
            goodRows(ii) = sum(isnan(alWaves(ii,:))) ~= size(alWaves,2);
        end
        [coeffs,scores] = pca(alWaves(:,logical(goodCols)),'NumComponents',2);
        kDists = nan(maxK);
        for ik = 1:maxK,
            sumD = nan(ik,kReps);
            for ir = 1:kReps,
                [idx(:,ir),~,sumD(:,ir)]= kmeans(scores,ik);
%                 idx2(:,ir) = kmeans(alWaves(:,goodCols),ik,'Distance','correlation');
            end
            [ids(:,ik), modDists] = klModeKmeans(idx,sumD);
            kDists(1:ik,ik) = modDists;
%             ids2(:,ik) = klModeKmeans(idx2);
        end
        
%         [kHat,kGap,kErr,kStruct] = klGetGapv7(alWaves(:,goodCols),ids2,'reclust','kmeans','randreps',10,'-r','rand','randval','pca','-t','corr','ktype','correlation');
        if forceK,
            kHat = forceKVal;
        else
            [kHat,kGap,kErr,kStruct] = klGetGapv8(scores,ids,'reclust','kmeans','randreps',10,'-r','rand','randval','unif','-t','euc','ktype','sqeuclidean','-d',kDists,'perc',1,'nsig',1);
        end
        for ik = 1:kHat,
            tstCorr(ik) = corr(nanmean(alWaves(ids(:,kHat) == ik,goodCols),1)',signalSeed(goodCols)');
            grpMeans(ik,:) = nanmean(alWaves(ids(:,kHat) == ik,:),1);
        end
        outIDs = ids;
        outGrp = find(tstCorr == max(tstCorr));
        isAP = ids(:,kHat) == outGrp;
        isNoise = ~isAP;
        
    case {'classify'},
        [coeffs,scores,~,~,expl] = pca(alWaves(:,goodCols),'NumComponents',3);
        grp=classify(scores,[scores(sortInd(1:nTake),:);scores(sortInd((end-nTake+1):end),:)],[ones(nTake,1);ones(nTake,1).*2]);
        outIDs = grp;
        isAP = grp == 2;
        isNoise = grp == 1;
        outGrp = 2;
end

isAP = logical(isAP); isNoise = logical(isNoise);


apMean = nanmean(alWaves(isAP,:),1);
% while size(apMean,2) < length(smoothTimes),
%     apMean(end+1) = nan;
% end
apTimes = alTimes;

% % Try plotting to visualize
if visualize
%     scatter(wvCorrs(:,1),wvCorrs(:,2),[],'b');
%     hold on;
%     plot(0:.01:1,0:.01:1,'r');
%     
%     xlabel('Signal correlation');
%     ylabel('Noise correlation');
%     
%     keyboard
    figure();
    plot(alTimes,alWaves(~isAP,:),'color',[.2 .2 .8]);
    hold on;
    plot(alTimes,alWaves(isAP,:),'color',[.8 .2 .2]);
    
    keyboard
end
    
    
% end
% 
