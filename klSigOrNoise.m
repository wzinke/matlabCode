% % klAuditQualitiesv1

%% Get constants
xlFile = './klDataBookKeeping_mg.xlsx';
upsamp = 1;
visualize = 1;
distType = 'corr';
wvTimes = ((1:32)-9).*25;
outFold = './Plots/allSorts';
xlWrite = 0;
online = 1;
doLater = 1;

global excelNum excelAll
load xlCols
% Edit the following line later for easier addition of waves instead of
% all-out recalculating...
sortWaves = [];
startTic = tic;
colors = [.8 .2 .2; .2 .8 .2; .2 .2 .8; .8 .2 .8; .2 .8 .8; .5 .5 .5];

load allSorts.mat

monk = fieldnames(allSorts);


for im = 2:length(monk),
    myMonk = monk{im};
    [excelNum,~,excelAll] = xlsread(xlFile,myMonk);
    
    kCol = find(strcmpi(excelAll(4,:),'kGuess'),1);
    goodCol = find(strcmpi(excelAll(4,:),'good'),1);
    sigCol = find(strcmpi(excelAll(4,:),'isSig'),1);
    
    % Get good (i.e., sorted) channels with k=1
    isVect = find(excelNum(:,kCol) == 1 & excelNum(:,goodCol) == 1);
    
    sigVect = nan(length(isVect),1);
    for ii = 1:length(isVect),
        if isnan(sigVect(ii)),
            is = isVect(ii);

            [path,file] = klRowToFile(is,'-m',myMonk,'-r',0,'-w',1);
            load([path,file{1}]);
            
            randInd = randperm(size(wave.waves,1),min([size(wave.waves,1),10000]));
            subWaves = wave.waves(randInd,:);
            
            % Display figures
            fig1 = figure('position',[66 499 560 420]);
            plot(wvTimes,subWaves);
            set(gca,'XLim',[-300,600]);

            [coeffs,scores] = pca(subWaves);

            fig2 = figure('position',[66 28 560 420]);
            scatter(scores(:,1),scores(:,2),[],'k');

            sig = 0;
            fprintf('Showing %s: %d of %d\n',monk{im},ii,length(isVect));
            fprintf('\tIf noise, hit F5. If signal, type "sig=1"\n');
            keyboard

            sigVect(ii) = sig;
            close all
        end
    end
    
    preSig = excelNum(:,sigCol);
    preSig(excelNum(:,kCol) > 1) = 1;
    preSig(isVect) = sigVect;
    
    xlswrite(xlFile,preSig(5:end,:),myMonk,sprintf('%s5',num2abc(sigCol)));
end
