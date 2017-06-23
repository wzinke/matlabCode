function [waveClustIDs, endk, clustTimes, clustSDF] = klWaveAgglomv3(varargin)

close all;

startTic = tic;

% Set user-defined variables
k           = 2:10;
nMembers    = 10;
waveType    = 'corr';
factType    = 'euc';
compType    = 'all';
printProg   = 1;

refine      = 1;
refType     = 'same';
kType       = 'correlation'; % Use "correlation" or "sqeuclidean"
xlWrite     = 0;
plotClusts  = 1; if nargin > 0, plotClusts = 0; end
plotDend    = 0;
plotFacts   = 0;
plotGroups  = 0; if nargout > 2, plotGroups = 1; end
plotIndiv   = 0;
smoothWaves = 1;

reextract   = 0;
reClust     = 1;
reGroup     = 1;
selectK     = 1;

areas = {'FEF'};
cutFacts = {'Pos','Amp','Width','Range'};
alignEvents = {'StimOn','SaccEnd'};

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'-k'}
            k = varargin{varStrInd(iv)+1};
        case {'-r','reextract'}
            reextract = varargin{varStrInd(iv)+1};
        case {'-n','nmembs'}
            nMembers = varargin{varStrInd(iv)+1};
        case {'-s','smooth'}
            smoothWaves = varargin{varStrInd(iv)+1};
        case {'-a','area'}
            areas       = varargin{varStrInd(iv)+1};
        case {'select'}
            selectK     = varargin{varStrInd(iv)+1};
    end
end

% Set constants
monk = {'Gauss','Helmholtz'};
xlFile = 'klDataBookKeeping_mg.xlsx';
types  = {'vis','mov','vismov','none'};

wvSampFreq = 40000;
wvTimeInc  = 1/wvSampFreq;
wvTimes    = ((1:32).*25)-(9*25); 
splTimes   = ((1:.1:32).*25)-(9*25);

areaStr = '';
for ia = 1:length(areas),
    areaStr = cat(2,areaStr,areas{ia});
end

close all
global excelNum excelAll wvClusts catFacts factNames keepFacts rawWaves allWaves normWaves wvTimes monkID catCrit

% Load column map
makeColLookup
load xlCols

% Initialize some variables
loadTic = tic;
fprintf('Extracting waveforms and measured characteristics...');
if isempty(monkID), reextract = 1; end
if reextract || ~exist('allWaves','var') || isempty(allWaves),
    allWaves = []; monkID = []; catFacts = []; catCrit = [];
    for im = 1:length(monk),
        [excelNum,~,excelAll] = xlsread(xlFile,monk{im});

        if ismember(areas,'All'),
             myCrit = logical([zeros(4,1);ones((size(excelAll,1)-4),1)]);
        else
            myCrit = ones(size(excelAll,1),1);
            for ia = 1:length(areas),
                myCrit = myCrit & strcmp(excelAll(:,col.area),areas{ia});
            end
        end
        critInds    = find(myCrit);

        % Get waveforms
        monkWaves = klGetWaveforms(critInds,'-m',monk{im});
        allWaves  = cat(1,allWaves,monkWaves);

        % Get column map for spike measures
        wdCol       = col.spkStart + (find(strcmpi(col.spkNames,'spkWidth'),1)-1);
        tfrCol      = col.spkStart + (find(strcmpi(col.spkNames,'TfR'),1)-1);
        ampCol      = col.spkStart + (find(strcmpi(col.spkNames,'amplitude'),1)-1);
        relAmpCol   = ampCol + 1;
        mnRateCol   = col.spkStart + (find(strcmpi(col.spkNames,'meanRate'),1)-1);
        fanoCol     = col.spkStart + (find(strcmpi(col.spkNames,'fano'),1)-1);
        cvCol       = col.spkStart + (find(strcmpi(col.spkNames,'cv'),1)-1);
        cv2Col      = col.spkStart + (find(strcmpi(col.spkNames,'cv2'),1)-1);
        lvCol       = col.spkStart + (find(strcmpi(col.spkNames,'lv'),1)-1);
        lvrCol      = col.spkStart + (find(strcmpi(col.spkNames,'lvr'),1)-1);
        mnISICol    = col.spkStart + (find(strcmpi(col.spkNames,'meanISI'),1)-1);

        % Get spike measures from excel file
        rawWd       = excelNum(myCrit,wdCol);
        rawTFR      = excelNum(myCrit,tfrCol);
        rawAmp      = excelNum(myCrit,ampCol);
        rawRange    = excelNum(myCrit,relAmpCol);
        rawRate     = excelNum(myCrit,mnRateCol);
        rawFano     = excelNum(myCrit,fanoCol);
        rawCV       = excelNum(myCrit,cvCol);
        rawCV2      = excelNum(myCrit,cv2Col);
        rawLV       = excelNum(myCrit,lvCol);
        rawLVR      = excelNum(myCrit,lvrCol);
        rawISI      = excelNum(myCrit,mnISICol);
        rawPos      = rawAmp < 0;
        rawAmp      = abs(rawAmp);

        factMat = [rawPos,rawWd,rawAmp,rawRange,rawRate,rawFano,rawCV,rawCV2,rawLV,rawISI];
        factNames = {'Pos','Width','Amp','Range','MnRate','Fano','CV','CV2','LV','ISI'};

        monkID                  = cat(1,monkID,ones(size(factMat,1),1).*im);
        catFacts                = cat(1,catFacts,factMat);
        catCrit                 = cat(1,catCrit,critInds);
    end
    rawWaves = allWaves;
end
fprintf('Waveforms extracted in %s\n',printTiming(loadTic));

if smoothWaves 
    if size(allWaves,2) == 32,
        for ii = 1:size(allWaves,1),
            splWaves(ii,:) = spline(1:32,allWaves(ii,:),1:.1:32);
        end
        allWaves = splWaves;
    end
    wvTimes = splTimes;    
end

% Align on derivative
[allWaves, wvTimes] = klTroughAlign(allWaves,wvTimes,0);

keepFacts = ~ismember(factNames,cutFacts);
cutFacts = catFacts(:,keepFacts);

factMean = nanmean(cutFacts,1);
factStd  = nanstd(cutFacts,[],1);
normFacts = (cutFacts - repmat(factMean,size(cutFacts,1),1))./repmat(factStd,size(cutFacts,1),1);
normWaves = allWaves - repmat(rawWaves(:,1),1,size(allWaves,2));

if reClust,
    [waveClustIDs, wvSimInd, wvClusts, wvLink] = klAgglomv4(normWaves,'-k',k,'-r',refine,'-p',printProg,'-t',waveType,'-n',nMembers,'refType',refType,'kType',kType);
%     [factClustIDs, ~, ~, factLink] = klAgglomv2(normFacts,'-k',k','-r',refine,'-p',printProg,'-t',factType,'-c',compType);
end

endk = k;
if selectK && size(waveClustIDs,2) > 1,
    tmp = figure();
    subplot(2,1,1);
    plot(k,wvSimInd); title('"F" statistic by k-clusters');
    subplot(2,1,2);
    plot(k(1:(end-1))+1,diff(wvSimInd)); title('Change in "F" statistic');
    hline(0);
    thisK = input('Select value of K to continue: ');
    waveClustIDs = waveClustIDs(:,k == thisK);
    endk = thisK;
    close(tmp);
end

fprintf('\n*** Clustering Completed in %s  ***\n',printTiming(startTic));

% Plot waveforms and mean waveforms for the clusters
if plotClusts,
    clustInds = unique(waveClustIDs);
    endLegStr = {};
    for ic = 1:length(clustInds),
        figure(ic);
        waveInds = waveClustIDs == clustInds(ic);
        clustMeans(ic,:) = nanmean(allWaves(waveInds,:),1);
        
        plot(wvTimes,allWaves(waveInds,:)');
        hold on;
        plot(wvTimes,clustMeans(ic,:),'k','linewidth',3);
        t=title(sprintf('Cluster %d (n=%d)',ic,sum(waveInds)));
        endLegStr{end+1} = sprintf('Clus. %d (n=%d)',ic,sum(waveInds));
    end

    figure(100);
    plot(wvTimes,clustMeans,'linewidth',3);
    set(gca,'XTick',0:200:800,'fontsize',16);
    xlabel('Time (us)','fontsize',20);
    ylabel('Relative Voltage','fontsize',20);
    title('Mean Within Cluster Waveforms','fontsize',24);
    legend(endLegStr,'fontsize',16);
    
    if plotDend
        figure(1000);
        dendrogram(wvLink,0);
    end
end

if plotFacts,
    clustInds = unique(waveClustIDs);
    factInds = unique(factClustIDs);
    for ic = 1:length(clustInds),
        for f = 1:length(factInds),
            numOverlap(ic,f) = sum(waveClustIDs == clustInds(ic) & factClustIDs == factInds(f));
            percClust(ic,f)  = numOverlap(ic,f)/sum(waveClustIDs == clustInds(ic));
            percFact(ic,f)   = numOverlap(ic,f)/sum(factClustIDs == factInds(f));
        end
    end
    figure(200);
    bar(percClust,'stacked');
    title('Percent Overlap of Wave/Factor Clusters','fontsize',24);
    xlabel('Waveform Cluster','fontsize',20);
    ylabel('Percent Wave Cluster','fontsize',20);
    set(gca,'fontsize',16,'YTick',0:.2:1,'YTickLabel',0:20:100);
    
end
    
if plotGroups,
    clustInds = unique(waveClustIDs);
    for ic = 1:length(clustInds),
        fprintf('Plotting Cluster %d (%d members)\n',ic,sum(waveClustIDs == clustInds(ic)));
        if sum(waveClustIDs == clustInds(ic)) > 100, continue; end
        clustSDF{ic} = {[],[]};
        clustTimes{ic} = {[],[]};
        for im = 1:length(unique(monkID)),
            lastMonk = find(monkID == (im-1),1,'last');
            if isempty(lastMonk), lastMonk = 0; end;
            thisGrpInds = catCrit(waveClustIDs == clustInds(ic) & monkID == im);
            if reGroup
                [grpSDF, grpTimes] = klPlotGrp(thisGrpInds,'-m',monk{im},'-p',0);
            end
            for ia = 1:length(grpSDF),
                % Find the min and max across set of times
                minTime = min([clustTimes{ic}{ia},grpTimes{1,ia}]);
                maxTime = max([clustTimes{ic}{ia},grpTimes{1,ia}]);
                newTimes = minTime:1:maxTime;
                
                % Make a matrix of nans the size of the end-goal of this
                % step
                newSDF = nan(size(clustSDF{ic}{ia},1)+size(grpSDF{ia},1),length(newTimes));
                % Put in the previously assigned cluster values
                newSDF(1:size(clustSDF{ic}{ia},1),ismember(newTimes,clustTimes{ic}{ia})) = clustSDF{ic}{ia};
                newSDF((size(clustSDF{ic}{ia},1)+1):size(newSDF,1),ismember(newTimes,grpTimes{1,ia})) = grpSDF{ia};
                
                % Store these values in clustSDF and clustTimes
                clustSDF{ic}{ia} = newSDF;
                clustTimes{ic}{ia} = newTimes;
                
               
            end
%             keyboard
        end
        for ia = 1:length(clustSDF{ic}),
            figure(ic+10*ia);
            if plotIndiv
                plot(clustTimes{ic}{ia},clustSDF{ic}{ia}); hold on;
                plot(clustTimes{ic}{ia},nanmean(clustSDF{ic}{ia},1),'k','linewidth',3);
            else
                pltMeanStd(clustTimes{ic}{ia},nanmean(clustSDF{ic}{ia},1),nanstd(clustSDF{ic}{ia},[],1)./sqrt(size(clustSDF{ic}{ia},1)),'k');
            end
            xlabel(sprintf('Time From %s',alignEvents{ia}),'fontsize',18);
            ylabel('Z-Scored SDF','fontsize',18); 
            title(sprintf('Cluster %d - Aligned on %s',ic,alignEvents{ia}),'fontsize',24);
        end
    end
end
    
%     waveClustIDs
%     [grpSDF
% 
% Output clusters to excel
if xlWrite
    lastCol  = find(strcmp(excelAll(4,:),'lastClust'));
    colHead  = ['corrClust',areaStr];
    clustCol = find(strcmp(excelAll(4,:),colHead));
    for im = 1:length(monk),
        clustsWrite     = nan(sum(monkID == im),1);
        clustsWrite(catCrit(monkID == im)) = waveClustIDs(monkID == im);
        xlswrite(xlFile,clustsWrite,monk{im},sprintf('%s1',num2abc(lastCol)));
        xlswrite(xlFile,{'lastClust'},monk{im},sprintf('%s4',num2abc(lastCol)));
        if ~isempty(clustCol),
            xlswrite(xlFile,clustsWrite,monk{im},sprintf('%s1',num2abc(clustCol)));
            xlswrite(xlFile,{colHead},monk{im},sprintf('%s4',num2abc(clustCol)));
        else
            warning('Unable to find appropriate column. Written to "last" only');
        end
    end
end