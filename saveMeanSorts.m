monk = {'Gauss','Helmholtz','Darwin'};
load allSorts.mat
task = 'MG';

for im = 1:length(monk),
    monkSort = allSorts.(monk{im});
    fprintf('Saving row 1 of %d...',length(monkSort.sortInfo));
    backStr = sprintf('1 of %d...',length(monkSort.sortInfo));
    for ir = 1:length(monkSort.sortInfo),
        if isempty(monkSort.sortInfo(ir).k), continue; end
        if mod(ir,20) == 0,
            for ib = 1:length(backStr), fprintf('\b'); end
            backStr = sprintf('%d of %d...',ir,length(monkSort.sortInfo));
            fprintf('%s',backStr);
        end
        [path,file] = klRowToFile(ir,'-w',1,'-m',monk{im},'-t',task);
        if isempty(file), continue; end;
        load([path,file{1}]);
        
        myInfo = monkSort.sortInfo(ir);
        myWaves = wave.waves(myInfo.waveInds,:);
        apWaves = myWaves(myInfo.groups == myInfo.apGroup,:);
        
        diffWaves = [nan(size(subWaves,1),5),diff(subWaves(:,5:20),[],2),nan(size(subWaves,1),length(21:32))];
        cutWaves=subWaves;
        cutWaves(abs(diffWaves) < 2 & abs(subWaves) > nanstd(nanmean(subWaves,1))*2.5) = nan;
        smoothWaves = nan(size(subWaves,1),length(1:.1:32));
        for ii = 1:size(subWaves,1),
            smoothWaves(ii,:) = spline(1:32,cutWaves(ii,:),1:.1:32);
        end
        smoothTimes = spline(1:32,times,1:.1:32);
        

        splWaves = nan(size(apWaves,1),length(1:.1:32));
        for iw = 1:size(apWaves,1),
            splWaves(iw,:) = spline(1:32,apWaves(iw,:),1:.1:32);
        end
        splTimes = spline(1:32,wvTimes,1:.1:32);
        
        [alWaves,meanSorts.(monk{im}).time{ir}] = klTroughAlignv4(splWaves,splTimes,0);
        meanSorts.(monk{im}).wave{ir} = nanmean(alWaves,1);
    end
    fprintf('\n');
end

save meanSorts.mat meanSorts