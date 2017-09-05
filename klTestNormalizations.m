%% KL plot Z-score problem

function klTestNormalizations()
% I made this a function so I could include some of mine below... you
% should be able to run this as if it were a script

% Load the file
pathPref = 'Y:/Users/Wolf/'; % WZ TEBA folder
path = [pathPref, 'ephys_db/Gauss/2014-12-02a/DSP/DSP10a/'];
file = '2014-12-02a_DSP10a_MG_a';
fStruct = load([path,file]);
spiketimes = fStruct.spiketimes;
Task = fStruct.Task;

% Get SDF (abbreviated klSpkRatev2 included below)
[sdf,sdfTimes] = klSpkRatev2(spiketimes);

% Get baseline mean/median and SD/MAD and do normalizations
blMean = nanmean(nanmean(sdf(:,sdfTimes > -200  & sdfTimes < 0),1),2);
blStd = nanstd(nanmean(sdf(:,sdfTimes > -200  & sdfTimes < 0),1),[],2);
zsdf = (sdf - blMean)./blStd;

blMedian = nanmedian(nanmean(sdf(:,sdfTimes > -200  & sdfTimes < 0),2),1);
blMAD = mad(nanmean(sdf(:,sdfTimes > -200  & sdfTimes < 0),2),1,1);
msdf = (sdf - blMedian)./blMAD;

% New figure
figure()
% Plot mean +/- sem of the SDF
subplot(2,2,1);
pltMeanStd(sdfTimes,nanmean(sdf,1),nanstd(sdf,1)./sqrt(size(sdf,1)),'k');
set(gca,'XLim',[-400,1000],'YLim',[0 40])
set(gca,'XLim',[-500,1000],'YLim',[0 40])
title('Mean +/- SEM');

% Plot median +/- MAD/sqrt(n) of the SDF
subplot(2,2,2);
pltMeanStd(sdfTimes,nanmedian(sdf,1),mad(sdf,1,1)./sqrt(size(sdf,1)),'k');
set(gca,'XLim',[-500,1000],'YLim',[0 40])
set(gca,'XLim',[-500,1000],'YLim',[0 20])
title('Median +/- MAD/sqrt(n)');

% Plot the Z-Score (x-mu)/sigma where mu and sigma are derived from
% baseline data
subplot(2,2,3);
pltMeanStd(sdfTimes,nanmean(zsdf,1),nanstd(zsdf,1)./sqrt(size(zsdf,1)),'k');
set(gca,'XLim',[-500,1000],'YLim',[-1 6]);
hline(2);
title('(SDF - mu)/sigma');

% Plot the non-parametric "Z-Score" (x-median)/mad where median and mad are
% derived from baseline data
subplot(2,2,4);
pltMeanStd(sdfTimes,nanmean(msdf,1),nanstd(msdf,1)./sqrt(size(msdf,1)),'k');
set(gca,'XLim',[-500,1000],'YLim',[-1 8]);
hline(2)
title('(SDF - median)/MAD');



    function [spkRate, tRange, groupedSpikes] = klSpkRatev2(spikes,varargin)               % V2 - spikes,bin,step,tStart,tEnd,type)
        % Set defaults
        bin = 200; step = 50; tStart = 0; tEnd = 10000; type = 'sdf'; splitByGroups = 0;
        quiet = 0;
        kType = 'psp';
        gaussSD = 25;

        minSpkTime = min(spikes(1:numel(spikes)));
        movedSpks = spikes - repmat((minSpkTime - 1),size(spikes));
        newMax = max(movedSpks(1:numel(movedSpks)));
        tRange = [floor(minSpkTime), (floor(minSpkTime) + ceil(newMax) - 1)];
        spkInds = [];
        for ic = 1:size(spikes,3)
            spkInds = cat(3,spkInds,zeros(size(spikes,1),ceil(newMax)));
        end
        for ic = 1:size(spikes,3)
            for ir = 1:size(spikes,1)
                if sum(~isnan(movedSpks(ir,:,ic))) > 0,
                    spkInds(ir,ceil(movedSpks(ir,~isnan(movedSpks(ir,:,ic)),ic)),ic) = 1;
                    spkInds(ir,find(spkInds(ir,:,ic) == 1,1,'last')+1:end,ic) = nan;
                    spkInds(ir,1:(find(spkInds(ir,:,ic) == 1,1)-1),ic) = nan;
                else
                    spkInds(ir,:,ic) = nan;
                end
            end
        end
        kern = klGetKern(kType,'width',gaussSD);
        kern = kern.*1000; % Convert from spk/ms to spk/s

        %spkDensMat = convn(spkInds,kernGauss,'same');   % cuts off convolution tails
        for ic = 1:size(spikes,3)
            if ~quiet
                fprintf('\tConvolving channel %d of %d...\n',ic,size(spikes,3));
            end
            spkDensMat(:,:,ic) = conv2(spkInds(:,:,ic),kern,'same');   % cuts off convolution tails
        end
        spkRate = spkDensMat;
        %         for ic = 1:size(spikes,3),
        %             for ir = 1:size(spikes,1),
        %                 spkRate(ir,find(spkInds(ir,:,ic) == 1,1,'last')+1:end,ic) = nan;
        %             end
        %         end
        %spkDensMean = mean(spkDensMat,1);
        %spkDensSEM = std(spkDensMat,1)./size(spkDensMat,1); % SEM = sd/sqrt(n)

        for ib = 1:ceil(size(spkDensMat,2)/step)
            binStart = (ib-1)*step+1;
            binEnd = binStart + bin;

            spkRateTemp(:,ib,:) = nanmean(spkDensMat(:,max([binStart,1]):min([binEnd,size(spkDensMat,2)])),2);
            spkThisBin(:,ib,:) = ceil(spkRateTemp(:,ib).*bin);
        end  
        tRange = tRange(1):tRange(2);

    end

    function [line shape] = pltMeanStd(xVals,means,stds,color)

        alpha = .5;
        lineWidth = 1;

        upMeans = means+stds;
        downMeans = means-stds;
        if sum(isnan(means)) == length(means) || sum(isnan(stds)) == length(stds),
            means(length(means):max([length(means),length(stds)])) = nan;
            stds(length(stds):max([length(means),length(stds)])) = nan;
        end

        if any(isnan(upMeans)) || any(isnan(downMeans))
            warning('NaN detected in SD values...\n\tFill may not work properly\n');
            upMeans(isnan(stds))    = [];
            downMeans(isnan(stds))  = [];
            xVals(isnan(stds))      = [];
            means(isnan(stds))      = [];
            %keyboard;
        end

        xValsFill = [xVals,xVals(end:-1:1)];
        stdFill = [upMeans,downMeans(end:-1:1)];


        shape = fill(xValsFill,stdFill,color,'FaceAlpha',alpha,'edgecolor',color','edgealpha',alpha); hold on;
        line = plot(xVals,means,'color',color,'linewidth',lineWidth);
    end


end