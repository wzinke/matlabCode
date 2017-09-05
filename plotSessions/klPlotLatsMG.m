% klPlotLatsMG
clearvars; close all;

%% Set constants
monk = {'Gauss','Helmholtz'};
xlFile = 'klDataBookKeeping_mg.xlsx';
types = {'vis','mov','vismov'};
typeColors = [.8 .2 .2; .2 .2 .8; .2 .8 .2];
areaList = {};
outDir = 'C:\Users\Kaleb\Documents\MATLAB\Plots\Latency\CellTypes';
vBins = 0:5:200;
mBins = -500:15:0;
maxVLat = inf;

%% Load excel file and orient to columns
for im = 1:2,
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    posCol          = find(strcmp(excelAll(4,:),'amplitude'));
    depthCol        = find(strcmp(excelAll(4,:),'depth (chan#)'));
    areaCol         = find(strcmp(excelAll(4,:),'area'));
    typeCol         = find(strcmp(excelAll(4,:),'typeAlt'));
    vLatCol         = find(strcmp(excelAll(4,:),'visRise'));
    mLatCol         = find(strcmp(excelAll(4,:),'movRise'));
    
    %% Get area list
    theseAreas  = unique(excelAll(5:end,areaCol));
    hasQuest    = strfind(theseAreas,'?');
    questInd    = find(~cellfun(@isempty,hasQuest));
    for iq = 1:length(questInd)
        theseAreas{questInd(iq)}(hasQuest{questInd(iq)}) = [];
    end
    areaList    = cat(1,areaList,theseAreas(~ismember(theseAreas,areaList)));
    
    %% Start area loop
    for ia = 1:length(areaList),
        typeLatsV = cell(1,length(types));
        typeLatsM = cell(1,length(types));
        thisArea = strcmp(excelAll(:,areaCol),areaList{ia});
        %% Loop through cell types (vis, mov, vismov)
        for it = 1:length(types),
            % Get this type latencies
            figure(ia); subplot(1,2,im); hold on;
            thisType = strcmp(excelAll(:,typeCol),types{it});
            thisSetLatsV = excelNum(thisArea & thisType,vLatCol);
            thisSetLatsV(thisSetLatsV > maxVLat) = nan;
            typeLatsV{it} = sort(thisSetLatsV(~isnan(thisSetLatsV)));
            if ~isempty(typeLatsV),
                plot(typeLatsV{it},(1:length(typeLatsV{it}))./length(typeLatsV{it}),'color',typeColors(it,:),'linewidth',3);
            end
            vHist{ia}(it,:) = hist(typeLatsV{it},vBins);
            vLegStr{ia}{it} = [upper(types{it}(1)),lower(types{it}(2:end)),' - n=',num2str(length(typeLatsV{it}))];
            
            figure(ia+10); subplot(1,2,im); hold on;
            thisSetLatsM = excelNum(thisArea & thisType,mLatCol);
            typeLatsM{it} = sort(thisSetLatsM(~isnan(thisSetLatsM)));
            if ~isempty(typeLatsM),
                plot(typeLatsM{it},(1:length(typeLatsM{it}))./length(typeLatsM{it}),'color',typeColors(it,:),'linewidth',3);
            end
            mHist{ia}(it,:) = hist(typeLatsM{it},mBins);
            mLegStr{ia}{it} = [upper(types{it}(1)),lower(types{it}(2:end)),' - n=',num2str(length(typeLatsM{it}))];
            
        end
        %% Pretty up the CDFs
        figure(ia);
        t(ia) = title(['Monkey ',upper(monk{im}(1))]);
        xLab = xlabel('Visual Latency');
        l = legend(vLegStr{ia});
        if im == 1, 
            yLab = ylabel('Cumulative Density');
        else
            st=suptitle(sprintf('Visual Latency CDFs: Area %s',areaList{ia}));
        end
        if im == 2,
            set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
            areaTitle = areaList{ia}; areaTitle(ismember(areaTitle,'/\')) = '-';
            saveas(gcf,sprintf('%s/%s-VisLatCDF.png',outDir,areaTitle));
        end
        
        figure(ia+10);
        t(ia) = title(['Monkey ',upper(monk{im}(1))]);
        xLab = xlabel('Pre-Motor Latency');
        l = legend(mLegStr{ia});
        if im == 1, 
            yLab = ylabel('Cumulative Density');
        else
            st=suptitle(sprintf('Pre-Motor Latency CDFs: Area %s',areaList{ia}));
        end
        if im == 2,
            set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
            areaTitle = areaList{ia}; areaTitle(ismember(areaTitle,'/\')) = '-';
            saveas(gcf,sprintf('%s/%s-MovLatCDF.png',outDir,areaTitle));
        end
        
        %% Make histograms
        figure(ia+100); subplot(2,1,im);
        vBar{ia} = bar(vBins,vHist{ia}');
        for it = 1:length(types),
            set(vBar{ia}(it),'facecolor',typeColors(it,:));
        end
%         clear vChild
        
        figure(ia+110); subplot(2,1,im);
        mBar{ia} = bar(mBins,mHist{ia}');
        for it = 1:length(types),
            set(mBar{ia}(it),'facecolor',typeColors(it,:));
        end
%         clear mChild
        
        % Pretty up V Bar
        figure(ia+100); subplot(2,1,im);
        t=title(['Monkey ',upper(monk{im}(1))]);
        yLab = ylabel('Count');
        l=legend(vLegStr{ia});
        if im == 2,
            xLab = xlabel('Visual Latency');
            st   = suptitle(sprintf('Visual Latency Histogram: Area %s',areaList{ia}));
            set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
            areaTitle = areaList{ia}; areaTitle(ismember(areaTitle,'/\')) = '-';
            saveas(gcf,sprintf('%s/%s-VisLatHist.png',outDir,areaTitle));
        end
        
        % Pretty up M Bar
        figure(ia+110); subplot(2,1,im);
        t=title(['Monkey ',upper(monk{im}(1))]);
        yLab = ylabel('Count');
        l=legend(mLegStr{ia});
        if im == 2,
            xLab = xlabel('Pre-Motor Latency');
            st   = suptitle(sprintf('Pre-Motor Latency Histogram: Area %s',areaList{ia}));
            set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
            areaTitle = areaList{ia}; areaTitle(ismember(areaTitle,'/\')) = '-';
            saveas(gcf,sprintf('%s/%s-MovLatHist.png',outDir,areaTitle));
        end
        
    end
end