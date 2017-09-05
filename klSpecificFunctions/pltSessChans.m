function pltSessChans(monk,sess)

% pltSessChans
close all

%% Set some user-defined variables
%monk = 'Helmholtz';
%sess = '2015-03-12a';
alignEvents = {'StimOn','SRT'};
saveSess = 1;
saveChans = 1;
saveDir = './Plots/sessionSummary';

%% Set some constants
dataFold        = 'Y:/Users/Wolf/ephys_db';
xlFile          = './cosDataBookKeeping.xlsx';
colors          = 'krgbcm';

%% Initialize some variables


[excelNum,excelText,excelAll] = xlsread(xlFile,monk);
excelNum = cat(1,nan(1,size(excelNum,2)),excelNum);
excelNum = cat(2,nan(size(excelNum,1),size(excelAll,2)-size(excelNum,2)),excelNum);

sessCol  = find(strcmp(excelAll(1,:),'Session'));
chanCol  = find(strcmp(excelAll(1,:),'Depth'));
areaCol  = find(strcmp(excelAll(1,:),'Area'));
posCol = find(strcmp(excelAll(1,:),'posMax'));

theseRows = find(strcmp(excelAll(:,sessCol),sess));

area = excelAll{theseRows(1),areaCol};
areaTitle = area; areaTitle(ismember(areaTitle,'/\')) = '-'; areaTitle(ismember(areaTitle,'?')) = '';
isPos    = excelNum(theseRows,posCol);

if saveChans
    [chanMeans, chanStds, chanTimes, chanEvents] = klPlotChansv3(theseRows,'-omsc',0,monk,1,1);
else
    [chanMeans, chanStds, chanTimes, chanEvents] = klPlotChansv3(theseRows,'-omc',0,monk,1);
end

close all;
figure(1); 
for ie = 1:2
    thisChanUsed    = zeros(1,24);
    yLims = nan(size(chanMeans,1),2);
    xLims = nan(size(chanMeans,1),2);
    axH = depthAxes('-fwx',1,.25,.3*(ie-1)+.15);
    set(axH,'box','on')
    colorLabel = [];
    for iu = 1:size(chanMeans,1)
        axes(axH(excelNum(theseRows(iu),chanCol)));
        thisChanUsed(excelNum(theseRows(iu),chanCol)) = thisChanUsed(excelNum(theseRows(iu),chanCol)) + 1;
        pltMeanStd(chanTimes{iu,ie},chanMeans{iu,ie},chanStds{iu,ie},colors(thisChanUsed(excelNum(theseRows(iu),chanCol))));
        yLims(iu,:) = get(gca,'YLim');
        xLims(iu,:) = get(gca,'XLim');
        if isPos(iu), colorLabel = cat(2,colorLabel,excelNum(theseRows(iu),chanCol)); end
    end
    for iy = 1:24
        axes(axH(iy));
        yl = ylabel(num2str(iy),'fontsize',12);
        if ismember(iy,colorLabel), set(yl,'color','r'); end
        set(gca,'YTickLabel','','XTickLabel','');
        set(axH(iy),'XLim',[min(xLims(:,1)),max(xLims(:,2))],'XTick',min(xLims(:,1)):500:max(xLims(:,2)))
        if iy == 24, set(axH(iy),'XTickLabel',min(xLims(:,1)):1000:max(xLims(:,2)),'fontsize',12); xlabel('Time (ms)','fontsize',14,'fontweight','bold'); end
        vl(iy) = vline(0);
    end
    set(vl,'color','g','linestyle','--','linewidth',2);
    yExt = [min(yLims(:,1)),max(yLims(:,2))];
    t = title(axH(1),sprintf('Aligned on %s',alignEvents{ie}),'fontsize',14,'fontweight','bold');
    if ie == 1,
        yLab = suplabel('Channel Number','y');
        set(yLab,'fontsize',14,'fontweight','bold');
    end
    %set(axH,'YLim',yExt);
end
st = suptitle(sprintf('Monkey %s: %s (%s)',upper(monk(1)),sess,area));
set(st,'fontsize',18,'fontweight','bold');

if saveSess
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s/%s-%s.png',saveDir,monk,areaTitle,sess));
end