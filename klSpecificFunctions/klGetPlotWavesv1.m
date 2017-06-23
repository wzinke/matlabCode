allWaves = []; allSNR = []; reconstID = []; nSpks = []; trRate = []; blRate = [];
for iFile = 1:length(fNames),
    fprintf('Getting waves from file %d (of %d)...\n',iFile,length(fNames));
    chanDir = dir(sprintf('%s/%s/Channel*',procDir,fNames{iFile}));
    for ic = 1:length(chanDir),
        unitDir = dir(sprintf('%s/%s/Channel%d/Unit*',procDir,fNames{iFile},ic));
        for iu = 1:length(unitDir),
            if exist(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',procDir,fNames{iFile},ic,iu),'file')
                try
                    load(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',procDir,fNames{iFile},ic,iu));
                    allWaves = cat(1,allWaves,nanmean(spikes.waves,1));
                    allSNR = cat(2,allSNR,spikes.qualVect(1));
                    nSpks = cat(2,nSpks,size(spikes.waves,1));
                    reconstID = cat(1,reconstID,[iFile,ic,iu]);
                    fields = fieldnames(spikes.rateStruct.trial);
                    for ii = 1:length(fields),
                        rateVect(1,ii) = spikes.rateStruct.trial.(fields{ii});
                        rateVect(2,ii) = spikes.rateStruct.baseline.(fields{ii});
                    end
                    trRate = cat(1,trRate,rateVect(1,:));
                    blRate = cat(1,blRate,rateVect(2,:));
                    clear spikes
                catch
                    continue
                end
            end
        end
    end
end

        
myInds = allSNR >= 1.25;
visualize = 0;
startTic = tic;
vMean = []; vStd = [];
mMean = []; mStd = [];
allMnWv = []; allStdWv = [];
reconstSub = reconstID(myInds,:);
for ii = 1:size(reconstSub,1),
    fprintf('Checking unit %d of %d...\n',ii,size(reconstSub,1));
    load(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',procDir,fNames{reconstSub(ii,1)},reconstSub(ii,2),reconstSub(ii,3)));
    load(sprintf('%s/%s/Behav.mat',procDir,fNames{reconstSub(ii,1)}));
    
    if visualize,
        figure();
        subplot(1,2,1);
        pltMeanStd(vTimes,nanmean(vSDF,1),nanstd(vSDF,1)./sqrt(size(vSDF,1)));
        set(gca,'XLim',[-300,700]);
        subplot(1,2,2);
        pltMeanStd(mTimes,nanmean(mSDF,1),nanstd(mSDF,1)./sqrt(size(mSDF,1)));
        set(gca,'XLim',[-700,300]);

        wvAx = axes('position',[.7 .7 .2 .2]);
        errorbar(1:32,nanmean(spikes.waves,1),nanstd(spikes.waves,1));

        suptitle(sprintf('%s-Ch%d-U%d: SNR=%.2f',fNames{reconstSub(ii,1)},reconstSub(ii,2),reconstSub(ii,3),spikes.qualVect(1)));
        pause
    end
    
    try
        [vSDF,vTimes] = klSpkRatev2(spikes.spiketimes(ismember(Task.TaskType,'MG'),:),'-q',1);
        [mSDF,mTimes] = klSpkRatev2(spikes.spiketimes(ismember(Task.TaskType,'MG'),:)-repmat(Task.SRT(ismember(Task.TaskType,'MG'))+Task.GoCue(ismember(Task.TaskType,'MG')),1,size(spikes.spiketimes,2)),'-q',1);
        vMean = cat(1,vMean,nanmean(vSDF(:,ismember(vTimes,-300:700)),1));
        mMean = cat(1,mMean,nanmean(mSDF(:,ismember(mTimes,-700:300)),1));
        vStd = cat(1,vStd,nanstd(vSDF(:,ismember(vTimes,-300:700)),[],1)./size(vSDF,1));
        mStd = cat(1,mStd,nanstd(mSDF(:,ismember(mTimes,-700:300)),[],1)./size(mSDF,1));
        allMnWv = cat(1,allMnWv,nanmean(spikes.waves,1));
        allStdWv = cat(1,allStdWv,nanstd(spikes.waves,[],1));
    catch
        vMean = cat(1,vMean,nan(1,size(vMean,2)));
        mMean = cat(1,mMean,nan(1,size(mMean,2)));
        vStd = cat(1,vStd,nan(1,size(vStd,2)));
        mStd = cat(1,mStd,nan(1,size(mStd,2)));
        allMnWv = cat(1,allMnWv,nan(1,size(allMnWv,2)));
        allStdWv = cat(1,allStdWv,nan(1,size(allStdWv,2)));
    end
    
    if size(vMean,1) ~= size(mMean,1),
        fprintf('***Problem at ii=%d...\n',ii);
    end
    
    
    close all
end

fprintf('Completed in %s\n',printTiming(startTic));

close all;
plotTic = tic;
for ii = 1:size(reconstSub,1),
    figure(ii);
    subplot(1,2,1);
    pltMeanStd(-300:700,vMean(ii,:),vStd(ii,:));
    set(gca,'XLim',[-300,700]);
    subplot(1,2,2);
    pltMeanStd(-700:300,mMean(ii,:),mStd(ii,:));
    set(gca,'XLim',[-700,300]);

    wvAx = axes('position',[.7 .7 .2 .2]);
    errorbar(1:32,allMnWv(ii,:),allStdWv(ii,:));

    suptitle(sprintf('%s-Ch%d-U%d: SNR=%.2f',fNames{reconstSub(ii,1)},reconstSub(ii,2),reconstSub(ii,3),spikes.qualVect(1)));
%     pause
    saveas(gcf,sprintf('./tdtResponses/%s-Ch%d-U%d.png',fNames{reconstSub(ii,1)},reconstSub(ii,2),reconstSub(ii,3)));
    close all
end
fprintf('Plots saved in %s\n',printTiming(plotTic));