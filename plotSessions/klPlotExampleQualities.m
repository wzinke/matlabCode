close all; warning off;

whichQuals = .25:.05:1;
qualCol = 2;

monk = {'Gauss','Helmholtz'};
wvTimes = ((1:32)-9).*25;
xlFile = './klDataBookKeeping_mg.xlsx';
nPerQual = 5;

makeColLookup
load xlCols

allIso = []; rowID = []; monkID = [];
for im = 2:length(monk),
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    allIso = cat(1,allIso,excelNum(5:end,col.SNR:col.SNR+3));
    rowID = cat(1,rowID,(5:size(excelNum,1))');
    monkID = cat(1,monkID,ones(size(excelNum,1)-4,1).*im);
end

for iq = 1:length(whichQuals),
    % Get closest observation to whichQuals
    dist = abs(allIso(:,qualCol) - whichQuals(iq));
    [sortDist,sortInd] = sort(dist);
    myRows = sortInd(1:nPerQual);
    
    for ir = 1:nPerQual,
        fprintf('Plotting #%d closest to Iso=%.2f\n',ir,whichQuals(iq));
        myRow = myRows(ir);
        % Load in the file
        [path,file] = klRowToFile(rowID(myRow),'-m',monk{monkID(myRow)},'-r',1,'-w',1);
        load([path,file{1}]);

        % Sort APs
%         diffWaves = [nan(size(wave.waves,1),5),diff(wave.waves(:,5:20),[],2),nan(size(wave.waves,1),length(21:32))];
%         cutWaves=wave.waves;
%         cutWaves(abs(diffWaves) < 2 & abs(wave.waves) > nanstd(nanmean(wave.waves,1))*2.5) = nan;
%         smoothWaves = nan(size(wave.waves,1),size(1:.1:32,2));
%         for ii = 1:size(wave.waves,1),
%             smoothWaves(ii,:) = spline(1:32,cutWaves(ii,:),1:.1:32);
%         end
%         smoothTimes = spline(1:32,wvTimes,1:.1:32);
%         alignWind = 10;
% %         smoothWaves = wave.waves;
% %         smoothTimes = wvTimes;
% %         alignWind = 2;
%         [alWaves, alTimes] = klTroughAlignv2(smoothWaves,smoothTimes,0,'-w',alignWind);

        [isAP,isNoise, sortOut, apMean, apTimes] = klSortSpikesv3(wave.waves,'-u',1,'-a',1,'-t',wvTimes,'type','classify');


        figure((iq*10)+ir);
        plot(apTimes,sortOut(~isAP,:),'b');
        hold on;
        plot(apTimes,sortOut(isAP,:),'r');
%         legend('Noise','Signal');
        xlabel('Time (us)'); ylabel('Relative Voltage');
        title(sprintf('%s Row %d - Iso/Miss/FA = %.2f/%.2f/%.2f',monk{monkID(myRow)}(1),rowID(myRow),allIso(myRow,2),allIso(myRow,3),allIso(myRow,4)));
        
        set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(gcf,sprintf('./Plots/sortQualities/%s%d-Iso%d.png',monk{monkID(myRow)}(1),rowID(myRow),round(allIso(myRow,2)*100)));
        close gcf;
        
        
    end
    
end