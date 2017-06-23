monk = {'Gauss','Helmholtz'};
xlFile = './klDataBookKeeping_mg.xlsx';
colors = [.8 .2 .2; .2 .8 .2; .2 .2 .8; .8 .2 .8; .2 .8 .8; .5 .5 .5];

for im = 1:length(monk),
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});

    isoMin = .9;

    rows = find(excelNum(:,col.SNR+1) >= isoMin);

    for ir = 1:length(rows),
        fprintf('Loading row %d (%d of %d)\n',rows(ir),ir,length(rows));
        [path,file] = klRowToFile(rows(ir),'-m',monk{im},'-r',1,'-w',1);
        load([path,file{1}]);
        [isAP,isNoise,grpMeans,apTimes,outGrp,alWaves,grpIDs] = klSortSpikesv4(wave.waves,'type','kmeans');
        figure(ir); hold on;
        uIDs = unique(grpIDs(:,size(grpMeans,1)));
        for i = 1:length(uIDs),
            plot(apTimes,alWaves(grpIDs(:,size(grpMeans,1)) == uIDs(i),:),'color',colors(i,:));
        end
        plot(apTimes,alWaves(isAP,:),'k'); hold on;
        plot(apTimes,grpMeans(outGrp,:),'--r','linewidth',2);
        xlabel('Time (us)');
        ylabel('Relative Voltage');

        title(sprintf('Monkey %s: Unit %d - Iso: %.2f, SNR:%.2f, k=%d',monk{im}(1),rows(ir),excelNum(rows(ir),col.SNR+1),excelNum(rows(ir),col.SNR),size(grpMeans,1)));
%         keyboard

        set(ir,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(ir,sprintf('C:/Users/Kaleb/Documents/MATLAB/Plots/testSorts160219/%s%d-Iso%d-Waves.png',monk{im}(1),rows(ir),round(excelNum(rows(ir),col.SNR+1)*100)));
        
        goodCols = zeros(1,size(alWaves,2));
        for ii = 1:size(alWaves,2),
            goodCols(ii) = sum(~isnan(alWaves(:,ii))) ==size(alWaves,1);
        end
        
        [coeffs,scores] = pca(alWaves(:,logical(goodCols)));
        figure(ir*10); hold on;
        uIDs = unique(grpIDs(:,size(grpMeans,1)));
        for i = 1:length(uIDs),
            scatter(scores(grpIDs(:,size(grpMeans,1))==uIDs(i),1),scores(grpIDs(:,size(grpMeans,1))==uIDs(i),2),[],colors(i,:));
        end
        scatter(scores(isAP,1),scores(isAP,2),[],'k');
        
        set(ir*10,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(ir*10,sprintf('C:/Users/Kaleb/Documents/MATLAB/Plots/testSorts160219/%s%d-Iso%d-PCA.png',monk{im}(1),rows(ir),round(excelNum(rows(ir),col.SNR+1)*100)));
        
        close all
    end
end