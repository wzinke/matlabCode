%clearvars; close all;
function compareTargetSelection(targType,breakOut)

close all; 

% compareCellTypes
monk        = {'Gauss','Helmholtz'};                     % Cell array of strings for monkey names to analyze. Even if it's one monkey, still keep it a cell
plotType    = 'sign';
doWZ        = 1;

%% Define some constants and initialize counters
xlFile      = 'cosDataBookKeeping_capture.xlsx';
areaList    = {};
types       = {'none','vis','mov','vismov'};
statTypes   = {'vTrans','vSust','preSacc','postSacc'};
tLocs       = [0, 90, 180, 270];
%targType    = 'target';
%breakOut    = 'type';
baseFold     = 'C:\Users\Kaleb\Documents\MATLAB\Plots\percUnitTypes';
if strcmpi(targType,'target'), targFold = 'TargetLocation'; elseif strcmpi(targType,'trial'), targFold = 'TrialType'; end
if strcmpi(breakOut,'type'), splitFold = 'ByType'; elseif strcmpi(breakOut,'sign'), splitFold = 'BySign'; end
outFold = sprintf('%s/%s/%s',baseFold,targFold,splitFold);

for im = 1:2
    [excelNum,excelText,excelAll] = xlsread(xlFile,monk{im});
    excelNum = cat(1,nan(size(excelAll,1)-size(excelNum,1),size(excelNum,2)),excelNum);
    excelNum = cat(2,nan(size(excelNum,1),size(excelAll,2)-size(excelNum,2)),excelNum);
    
    % Find the header row and define columns
    hRow = find(strcmp(excelAll(:,1),'Session'),1);
    
    posCol          = find(strcmp(excelAll(hRow,:),'posMax'));
    depthCol        = find(strcmp(excelAll(hRow,:),'Depth'));
    areaCol         = find(strcmp(excelAll(hRow,:),'Area'));
    typeCol         = find(strcmp(excelAll(hRow,:),'Type'));
    rfCol           = find(strcmp(excelAll(hRow,:),'RF'));
    allVT           = find(strcmp(excelAll(hRow,:),'vTrans'));
    
    vtTarget = allVT(1);
    vtType   = allVT(2:end);
    
    isPos    = excelNum(:,posCol) == 1;
    isNeg    = excelNum(:,posCol) == 0;
    
    whichRF  = excelNum(:,rfCol);
    isCompared = ~isnan(whichRF);
    
    % Get unique area names
    theseAreas  = unique(excelAll((hRow+1):end,areaCol));
    hasQuest    = strfind(theseAreas,'?');
    questInd    = find(~cellfun(@isempty,hasQuest));
    for iq = 1:length(questInd)
        theseAreas{questInd(iq)}(hasQuest{questInd(iq)}) = [];
    end
    areaList    = cat(1,areaList,theseAreas(~ismember(theseAreas,areaList)));
    
    allComps = nan(size(excelAll,1),4);
    switch lower(targType)
        case 'target'
            allComps = excelNum(:,vtTarget:vtTarget+3);
        case 'trial'
            %% Loop through all rows
            for ir = (hRow+1):size(excelAll,1),
                % Isolate the comparisons done in the neuron's RF
                whichSet = find(ismember(tLocs,excelNum(ir,rfCol)));
                theseComps = excelNum(ir,vtType(whichSet):(vtType(whichSet)+3));
                if ~isempty(theseComps)
                    allComps(ir,:) = theseComps;
                end
            end
    end
    
    switch breakOut
        case 'type'
            %% Loop through areas
            for ia = 1:length(areaList)
                thisArea = strcmp(excelAll(:,areaCol),areaList{ia});
                %% Loop through types
                for it = 1:length(types)
                    thisType = strcmp(excelAll(:,typeCol),types{it});
                    sumTargetLoc{ia}(im,1:4,it) = repmat(sum(thisArea & thisType & isCompared),1,4);
                    %% Loop through tests to compare
                    for ip = 1:4
                        numSelTargetLoc{ia}(im,ip,it) = sum(thisArea & thisType & isCompared & allComps(:,ip) < .05);
                    end
                end
                percSelTargetLoc{ia} = numSelTargetLoc{ia}./sumTargetLoc{ia};
            end
        case 'sign'
            for ia = 1:length(areaList)
                thisArea = strcmp(excelAll(:,areaCol),areaList{ia});
                for is = 0:1
                    thisSign = excelNum(:,posCol) == is;
                    sumTargetLoc{ia}(im,1:4,is+1) = repmat(sum(thisArea & thisSign & isCompared),1,4);
                    %% Loop through tests to compare
                    for ip = 1:4
                        numSelTargetLoc{ia}(im,ip,is+1) = sum(thisArea & thisSign & isCompared & allComps(:,ip) < .05);
                    end
                end
                percSelTargetLoc{ia} = numSelTargetLoc{ia}./sumTargetLoc{ia};
                if strcmpi(breakOut,'sign'),
                    figure(ia); subplot(1,2,im); hold on;
                    bar([percSelTargetLoc{ia}(im,:,1);percSelTargetLoc{ia}(im,:,2)]');
                    t = title(sprintf('Monkey %s',upper(monk{im}(1))),'fontsize',12,'fontweight','bold');
                    set(gca,'YLim',[0 1],'YTick',0:.25:1,'YTickLabel',0:25:100,'XLim',[0 5], 'XTick',1:4,'XTickLabel',statTypes);
                    nNeg = sum(thisArea & isNeg); nPos = sum(thisArea & isPos);
                    l = legend({sprintf('Neg Spks (n = %d)',nNeg),sprintf('Pos Spks (n = %d)',nPos)});
                    
                    if im == 2,
                        areaTitle = areaList{ia};
                        areaTitle(ismember(areaTitle,'/\')) = '-';
                        areaTitle(ismember(areaTitle,'?')) = '';
                        switch targType
                            case 'target'
                                st = suptitle(sprintf('%s - Target Location Selectivity',areaList{ia}));
                            case 'trial'
                                st = suptitle(sprintf('%s - Trial Type Selectivity',areaList{ia}));
                        end
                        set(st,'fontsize',14,'fontweight','bold');
                        xx = suplabel('Comparison Type','x');
                        yy = suplabel('Percent Significant','y');
                        set([xx,yy],'fontsize',12,'fontweight','bold');
                        
                        set(ia,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
                        saveas(ia,sprintf('%s/%s-%s.png',outFold,areaTitle,targFold));
                    end
                    
                end
            end
    end
    %keyboard
end

for ia = 1:length(areaList)
    figure(ia);
    areaTitle = areaList{ia};
    areaTitle(ismember(areaTitle,'/\')) = '-';
    areaTitle(ismember(areaTitle,'?')) = '';
    if strcmpi(breakOut,'type'),
        for it = 1:4
            clear x y
            subplot(2,2,it);
            bar(percSelTargetLoc{ia}(:,:,it)'); 
            t = title(types{it},'fontsize',14,'fontweight','bold');
            set(gca,'YLim',[0 1],'YTick',0:.25:1,'YTickLabel',0:25:100,'XLim',[0 5], 'XTick',1:4,'XTickLabel',statTypes);
            if it == 2,
                l = legend('Gauss','Helmholtz','location','best');
                set(l,'box','off');
            end
        end
        switch targType
            case 'target'
                st = suptitle(sprintf('%s - Target Location Selectivity',areaList{ia}));
            case 'trial'
                st = suptitle(sprintf('%s - Trial Type Selectivity',areaList{ia}));
        end
        set(st,'fontsize',16,'fontweight','bold');
        xx = suplabel('Comparison Type','x');
        yy = suplabel('Percent Significant','y');
        set([xx,yy],'fontsize',12,'fontweight','bold');

        set(ia,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(ia,sprintf('%s/%s-%s.png',outFold,areaTitle,targFold));
    end
end

%keyboard