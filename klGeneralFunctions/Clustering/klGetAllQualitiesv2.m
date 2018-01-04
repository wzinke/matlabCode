% klGetAllQualities

% monk = {'Gauss','Helmholtz'};
monk = {'Darwin'};
xlFile = './klDataBookKeeping_mg.xlsx';
upsamp = 1;
visualize = 1;
distType = 'corr';
wvTimes = ((1:32)-9).*25;
outFold = './Plots/allSorts';
xlWrite = 0;
minWaves = 250;

global excelNum excelAll
load xlCols
% Edit the following line later for easier addition of waves instead of
% all-out recalculating...
sortWaves = [];
startTic = tic;
colors = [.8 .2 .2; .2 .8 .2; .2 .2 .8; .8 .2 .8; .2 .8 .8; .5 .5 .5];

for im = 1:length(monk)
    monkTic = tic;
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    qualMat = nan(size(excelAll,1)-4,5);
    startRow = 1;
    for ir = 5:size(excelAll,1),
        if mod(ir,1) == 0, fprintf('Sorting wave (%d/%d)\n',ir,size(excelAll,1)); end
        [path,file] = klRowToFile(ir,'-m',monk{im},'-w',1);
        load([path,file{1}]);
		
		% Changed to klUnitIsolationv3 on 2/19/16 (v2 should still work as before)
        if size(wave.waves,1) >= minWaves,
            [snr,isoScore,fnScore,fpScore,outK, outInfo] = klUnitIsolationv3(wave.waves,'-u',upsamp,'-d',distType,'-v',visualize,'-t',wvTimes);
        else
            snr = nan;
            isoScore = nan;
            fnScore = nan;
            fpScore = nan;
            outK = 1;
            outInfo.k = 1;
            outInfo.apGroup = 1;
            outInfo.groups = ones(size(wave.waves,1),1);
            outInfo.waveInds = 1:size(wave.waves,1);
        end
        
%         qualMat(ir-4,:) = [snr,isoScore,fnScore,fpScore,outK];
        
        if visualize
            thisFig = figure('visible','off');
            for ik = 1:outInfo.k
                plot(outInfo.times,outInfo.alWaves(outInfo.groups == ik,:),'color',colors(ik,:));
                hold on;
            end
            plot(outInfo.times,nanmean(outInfo.alWaves(outInfo.groups == outInfo.apGroup,:),1),'k','linewidth',2);
            title(sprintf('%s%d: SNR=%.2f, Iso=%.3f, FN=%.3f, FP=%.3f',monk{im}(1),ir,snr,isoScore,fnScore,fpScore));
            set(gca,'XLim',[-300,600]);
            
            set(thisFig,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
            saveas(thisFig,sprintf('%s/%s%d-Iso%d-Waves.png',outFold,monk{im}(1),ir,round(isoScore*100)));
            close(thisFig);
            
            for ii = 1:size(outInfo.alWaves,2),
                goodCols(ii) = sum(~isnan(outInfo.alWaves(:,ii))) == size(outInfo.alWaves(:,ii),1);
            end
            [coeffs,scores] = pca(outInfo.alWaves(:,logical(goodCols)));
            
            thisFig = figure('visible','off');
            for ik = 1:outInfo.k
                scatter(scores(outInfo.groups == ik,1),scores(outInfo.groups == ik,2),[],colors(ik,:));
                hold on;
            end
            scatter(scores(outInfo.groups == outInfo.apGroup,1),scores(outInfo.groups == outInfo.apGroup,2),[],'k');
            title(sprintf('%s%d: SNR=%.2f, Iso=%.3f, FN=%.3f, FP=%.3f',monk{im}(1),ir,snr,isoScore,fnScore,fpScore));
%             set(gca,'XLim',[-300,600]);
            
            set(thisFig,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
            saveas(thisFig,sprintf('%s/%s%d-Iso%d-PCA.png',outFold,monk{im}(1),ir,round(isoScore*100)));
            close(thisFig);
        end
        
        sortQualities.(monk{im}).sortInfo(ir).k = outInfo.k;
%         sortQualities.(monk{im}).sortInfo(ir).times = outInfo.times;
%         sortQualities.(monk{im}).sortInfo(ir).waves = outInfo.alWaves;
        sortQualities.(monk{im}).sortInfo(ir).apGroup = outInfo.apGroup;
        sortQualities.(monk{im}).sortInfo(ir).groups = outInfo.groups;
        sortQualities.(monk{im}).sortInfo(ir).waveInds = outInfo.waveInds;
        
        sortQualities.(monk{im}).quals(ir,:) = [snr,isoScore,fnScore,fpScore,outK];
    
        if mod(ir,100)==0 || ir == size(excelAll,1),
            save(sprintf('%sQualities_%d-%d.mat',monk{im}(1),startRow,ir),'sortQualities','-v7.3');
            startRow = ir+1;
            clear sortQualities
        end
    end
%     sortQualities.(monk{im}).quals = qualMat;
    
    if xlWrite,
        xlswrite(xlFile,sortQualities.(monk{im}),monk{im},sprintf('%s5',num2abc(col.SNR)));
    end
    fprintf('\nMonkey %s sorted in %s\n',monk{im},printTiming(monkTic));
end


fprintf('Completed in %s\n',printTiming(startTic));