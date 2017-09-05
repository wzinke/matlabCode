% function klSortAllSpikes

monk = {'Gauss','Helmholtz'};
xlFile = './klDataBookKeeping_mg.xlsx';
upsamp = 1;
sortType = 'kmeans';
visualize = 0;
wvTimes = ((1:32).*25)-(9*25);

global excelNum excelAll

% Edit the following line later for easier addition of waves instead of
% all-out recalculating...
fprintf('**** Starting to sort all waveforms. Will save to sortedWaves.mat  *****\n');
sortWaves = [];
allTic = tic;
for im = 1:length(monk)
    monkTic = tic;
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    
    for ir = 5:size(excelAll,1),
        if mod(ir,10) == 0, fprintf('\tSorting wave (%d/%d)\n',ir,size(excelAll,1)); end
        [path,file] = klRowToFile(ir,'-m',monk{im},'-w',1);
        load([path,file{1}]);
        [isAP,isNoise,grpMeans,apTimes,outGrp,alWaves,grpIDs] = klSortSpikesv3(wave.waves,'type','kmeans');
        
        sortWaves.(monk{im}).(sprintf('s%s',excelAll{ir,2}(~ismember(excelAll{ir,2},'-')))).(excelAll{ir,5}).apWave = grpMeans(outGrp,:);
        sortWaves.(monk{im}).(sprintf('s%s',excelAll{ir,2}(~ismember(excelAll{ir,2},'-')))).(excelAll{ir,5}).apTimes = apTimes;
    end
    save(sprintf('sortedWaves_%s.mat',monk{im}),'sortWaves');
    fprintf('Monkey %s Completed in %s\n',monk{im},printTiming(monkTic));
end

save('sortedWaves.mat','sortWaves');
fprintf('All waveforms sorted and saved in %s\n',printTiming(allTic));