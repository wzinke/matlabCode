% klPlotLats
close all;

saveDir = './Plots/Latency';
monk   = {'Gauss','Helmholtz'};
myArea = 'All';
xlFile = './klDataBookKeeping_mg.xlsx';

for im = 1:2,
    close all;
    clearvars -except im monk myArea saveDir xlFile
    [excelNum,excelText,excelAll] = xlsread(xlFile,monk{im});
    excelNum = cat(1,nan(4,size(excelNum,2)),excelNum);
    excelNum = cat(2,nan(size(excelNum,1),2),excelNum);

    hRow = find(strcmp(excelAll(:,2),'Name'),1);

    visLatCol   = find(strcmp(excelAll(hRow,:),'visLat'));
    movLatCol   = find(strcmp(excelAll(hRow,:),'movLat'));
    areaCol     = find(strcmp(excelAll(hRow,:),'area'));
    typeCol     = find(strcmp(excelAll(hRow,:),'type'));
    ampCol      = find(strcmp(excelAll(hRow,:),'amplitude'));

    if strcmpi(myArea,'all'),
        isMyArea = ones(size(excelAll,1),1);
    else
        isMyArea    = strcmpi(excelAll(:,areaCol),myArea);
    end
    isPos       = excelNum(:,ampCol) < 1;

    isVis       = strcmpi(excelAll(:,typeCol),'vis');
    isMov       = strcmpi(excelAll(:,typeCol),'mov');
    isVM        = strcmpi(excelAll(:,typeCol),'vismov');

    visLatsV    = excelNum(isMyArea & isVis,visLatCol);
    movLatsV    = excelNum(isMyArea & isMov,visLatCol);
    vmLatsV     = excelNum(isMyArea & isVM, visLatCol);
    visLatsM    = excelNum(isMyArea & isVis,movLatCol);
    movLatsM    = excelNum(isMyArea & isMov,movLatCol);
    vmLatsM     = excelNum(isMyArea & isVM, movLatCol);

    visLatsVPos    = excelNum(isMyArea & isVis & isPos,visLatCol);
    movLatsVPos    = excelNum(isMyArea & isMov & isPos,visLatCol);
    vmLatsVPos     = excelNum(isMyArea & isVM & isPos, visLatCol);
    visLatsMPos    = excelNum(isMyArea & isVis & isPos,movLatCol);
    movLatsMPos    = excelNum(isMyArea & isMov & isPos,movLatCol);
    vmLatsMPos     = excelNum(isMyArea & isVM & isPos, movLatCol);

    visLatsVNeg    = excelNum(isMyArea & isVis & ~isPos,visLatCol);
    movLatsVNeg    = excelNum(isMyArea & isMov & ~isPos,visLatCol);
    vmLatsVNeg     = excelNum(isMyArea & isVM & ~isPos, visLatCol);
    visLatsMNeg    = excelNum(isMyArea & isVis & ~isPos,movLatCol);
    movLatsMNeg    = excelNum(isMyArea & isMov & ~isPos,movLatCol);
    vmLatsMNeg     = excelNum(isMyArea & isVM & ~isPos, movLatCol);

    allPosV     = excelNum(isMyArea & isPos,visLatCol);
    allPosM     = excelNum(isMyArea & isPos,movLatCol);
    allNegV     = excelNum(isMyArea & ~isPos,visLatCol);
    allNegM     = excelNum(isMyArea & ~isPos,movLatCol);
    allPosV     = sort(allPosV(isfinite(allPosV)));
    allPosM     = sort(allPosM(isfinite(allPosM)));
    allNegV     = sort(allNegV(isfinite(allNegV)));
    allNegM     = sort(allNegM(isfinite(allNegM)));


    vBins = -100:5:600;
    mBins = -500:5:200;

    vMat(:,1) = hist(visLatsV,vBins);
    vMat(:,2) = hist(movLatsV,vBins);
    vMat(:,3) = hist(vmLatsV,vBins);

    vMatPos(:,1) = hist(visLatsVPos,vBins);
    vMatPos(:,2) = hist(movLatsVPos,vBins);
    vMatPos(:,3) = hist(vmLatsVPos,vBins);
    vMatPos(:,4) = hist(allPosV,vBins);

    vMatNeg(:,1) = hist(visLatsVNeg,vBins);
    vMatNeg(:,2) = hist(movLatsVNeg,vBins);
    vMatNeg(:,3) = hist(vmLatsVNeg,vBins);
    vMatNeg(:,4) = hist(allNegV,vBins);

    mMatPos(:,1) = hist(visLatsMPos,mBins);
    mMatPos(:,2) = hist(movLatsMPos,mBins);
    mMatPos(:,3) = hist(vmLatsMPos,mBins);
    mMatPos(:,4) = hist(allPosM,mBins);

    mMatNeg(:,1) = hist(visLatsMNeg,mBins);
    mMatNeg(:,2) = hist(movLatsMNeg,mBins);
    mMatNeg(:,3) = hist(vmLatsMNeg,mBins);
    mMatNeg(:,4) = hist(allNegM,mBins);

    mMat(:,1) = hist(visLatsM,mBins);
    mMat(:,2) = hist(movLatsM,mBins);
    mMat(:,3) = hist(vmLatsM,mBins);

    visDiff  = visLatsV-visLatsM; visDiff = visDiff(isfinite(visDiff));
    movDiff  = movLatsV-movLatsM; movDiff = movDiff(isfinite(movDiff));
    vmDiff   = vmLatsV-vmLatsM;   vmDiff  = vmDiff(isfinite(vmDiff));

    visLatsV = sort(visLatsV(isfinite(visLatsV)));
    movLatsV = sort(movLatsV(isfinite(movLatsV)));
    vmLatsV  = sort(vmLatsV(isfinite(vmLatsV)));

    visLatsM = sort(visLatsM(isfinite(visLatsM)));
    movLatsM = sort(movLatsM(isfinite(movLatsM)));
    vmLatsM  = sort(vmLatsM(isfinite(vmLatsM)));

    visLatsVPos = sort(visLatsVPos(isfinite(visLatsVPos)));
    movLatsVPos = sort(movLatsVPos(isfinite(movLatsVPos)));
    vmLatsVPos  = sort(vmLatsVPos(isfinite(vmLatsVPos)));

    visLatsMPos = sort(visLatsMPos(isfinite(visLatsMPos)));
    movLatsMPos = sort(movLatsMPos(isfinite(movLatsMPos)));
    vmLatsMPos  = sort(vmLatsMPos(isfinite(vmLatsMPos)));

    visLatsVNeg = sort(visLatsVNeg(isfinite(visLatsVNeg)));
    movLatsVNeg = sort(movLatsVNeg(isfinite(movLatsVNeg)));
    vmLatsVNeg  = sort(vmLatsVNeg(isfinite(vmLatsVNeg)));

    visLatsMNeg = sort(visLatsMNeg(isfinite(visLatsMNeg)));
    movLatsMNeg = sort(movLatsMNeg(isfinite(movLatsMNeg)));
    vmLatsMNeg  = sort(vmLatsMNeg(isfinite(vmLatsMNeg)));

    legStr = {sprintf('Vis: n=%d',length(visLatsV)),sprintf('Mov: n=%d',length(movLatsV)),sprintf('VisMov: n=%d',length(vmLatsV))};

    figure();
    bar(vBins,vMat);
    legend(legStr);
    t=title([monk{im},': Visual Latencies']);
    xlabel('Time After Stim (ms)');
    ylabel('# of cells');
    sortHistV = sort(vMat(:),'ascend');
    cutMax = sortHistV(find(diff(sortHistV) < 20,1,'last')+1);
    set(gca,'YLim',[0,cutMax+5]);
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s-CellTypeVisHist.png',saveDir,monk{im}));
    
    figure();
    bar(mBins,mMat);
    legend(legStr);
    t=title([monk{im},': Pre-Saccade Latencies']);
    xlabel('Time to Saccade (ms)');
    ylabel('# of cells');
    % Get reasonable ymax... 
    sortHistM = sort(mMat(:),'ascend');
    cutMax = sortHistM(find(diff(sortHistM) < 20,1,'last')+1);
    set(gca,'YLim',[0,cutMax+5]);
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s-CellTypeMovHist.png',saveDir,monk{im}));
    
    figure(); hold on;
    % klCDF(visLatsV,'color','k','plot',1);
    % klCDF(movLatsV,'color','r','plot',1);
    % klCDF(vmLatsV,'color','b','plot',1);
    plot(sort(visLatsV,'ascend'),(1:length(visLatsV))./length(visLatsV),'k');
    plot(sort(movLatsV,'ascend'),(1:length(movLatsV))./length(movLatsV),'r');
    plot(sort(vmLatsV,'ascend'),(1:length(vmLatsV))./length(vmLatsV),'b');
    legend(legStr,'location','nw');
    t=title([monk{im},': Visual Latencies']);
    set(gca,'YLim',[0,1],'YTick',0:.2:1,'YTickLabel',0:20:100);
    a(1) = xlabel('Visual Latency'); a(2) = ylabel('Cumulative Density (%)');
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s-CellTypeVisLat.png',saveDir,monk{im}));
    
    figure(); hold on;
    % klCDF(visLatsM,'color','k','plot',1);
    % klCDF(movLatsM,'color','r','plot',1);
    % klCDF(vmLatsM,'color','b','plot',1);
    plot(sort(visLatsM,'ascend'),(1:length(visLatsM))./length(visLatsM),'k');
    plot(sort(movLatsM,'ascend'),(1:length(movLatsM))./length(movLatsM),'r');
    plot(sort(vmLatsM,'ascend'),(1:length(vmLatsM))./length(vmLatsM),'b');
    legend(legStr,'location','nw');
    t=title([monk{im},': Pre-Saccade Latencies']);
    set(gca,'YLim',[0,1],'YTick',0:.2:1,'YTickLabel',0:20:100);
    a(1) = xlabel('Visual Latency'); a(2) = ylabel('Cumulative Density (%)');
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s-CellTypeMovLat.png',saveDir,monk{im}));
    
    figure()
    hold on;
    plot(sort(visDiff),(1:length(visDiff))./length(visDiff),'k');
    plot(sort(movDiff),(1:length(movDiff))./length(movDiff),'r');
    plot(sort(vmDiff),(1:length(vmDiff))./length(vmDiff),'b');
    legend(legStr,'location','nw');
    t=title([monk{im},': VisLat - MovLat']);
    set(gca,'YLim',[0,1],'YTick',0:.2:1,'YTickLabel',0:20:100);
    a(1) = xlabel('VisLatency - MovLatency'); a(2) = ylabel('Cumulative Density (%)');
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s-CellTypeLatDiff.png',saveDir,monk{im}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure();
    subplot(2,2,1); hold on;
    plot(visLatsVPos,(1:length(visLatsVPos))./length(visLatsVPos),'r');
    plot(visLatsVNeg,(1:length(visLatsVNeg))./length(visLatsVNeg),'k');
    t(1) = title('Vis Cells');
    legend(['Pos: n=',num2str(length(visLatsVPos))],['Neg: n=',num2str(length(visLatsVNeg))]);
    set(gca,'YLim',[0 1],'YTick',0:.5:1,'YTickLabel',0:50:100);

    subplot(2,2,2); hold on;
    plot(movLatsVPos,(1:length(movLatsVPos))./length(movLatsVPos),'r');
    plot(movLatsVNeg,(1:length(movLatsVNeg))./length(movLatsVNeg),'k');
    t(2) = title('Mov Cells');
    legend(['Pos: n=',num2str(length(movLatsVPos))],['Neg: n=',num2str(length(movLatsVNeg))]);
    set(gca,'YLim',[0 1],'YTick',0:.5:1,'YTickLabel',0:50:100);

    subplot(2,2,3); hold on;
    plot(vmLatsVPos,(1:length(vmLatsVPos))./length(vmLatsVPos),'r');
    plot(vmLatsVNeg,(1:length(vmLatsVNeg))./length(vmLatsVNeg),'k');
    t(3) = title('VisMov Cells');
    legend(['Pos: n=',num2str(length(vmLatsVPos))],['Neg: n=',num2str(length(vmLatsVNeg))]);
    set(gca,'YLim',[0 1],'YTick',0:.5:1,'YTickLabel',0:50:100);

    subplot(2,2,4); hold on;
    plot(allPosV,(1:length(allPosV))./length(allPosV),'r');
    plot(allNegV,(1:length(allNegV))./length(allNegV),'k');
    t(4) = title('All Cells');
    set(gca,'YLim',[0 1],'YTick',0:.5:1,'YTickLabel',0:50:100);
    legend(['Pos: n=',num2str(length(allPosV))],['Neg: n=',num2str(length(allNegV))]);
    st      = suptitle(['Visual Latencies - Monkey ',upper(monk{im}(1))]);
    sl(1)   = suplabel('Visual Latency','x');
    sl(2)   = suplabel('Cumulative Density (%)');
    set(st,'fontsize',16); set(sl,'fontsize',14);
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s-PosSpksVisLat.png',saveDir,monk{im}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure();
    subplot(2,2,1); hold on;
    plot(visLatsMPos,(1:length(visLatsMPos))./length(visLatsMPos),'r');
    plot(visLatsMNeg,(1:length(visLatsMNeg))./length(visLatsMNeg),'k');
    t(1) = title('Vis Cells');
    legend(['Pos: n=',num2str(length(visLatsMPos))],['Neg: n=',num2str(length(visLatsMNeg))]);
    set(gca,'YLim',[0 1],'YTick',0:.5:1,'YTickLabel',0:50:100);

    subplot(2,2,2); hold on;
    plot(movLatsMPos,(1:length(movLatsMPos))./length(movLatsMPos),'r');
    plot(movLatsMNeg,(1:length(movLatsMNeg))./length(movLatsMNeg),'k');
    t(2) = title('Mov Cells');
    legend(['Pos: n=',num2str(length(movLatsMPos))],['Neg: n=',num2str(length(movLatsMNeg))]);
    set(gca,'YLim',[0 1],'YTick',0:.5:1,'YTickLabel',0:50:100);

    subplot(2,2,3); hold on;
    plot(vmLatsMPos,(1:length(vmLatsMPos))./length(vmLatsMPos),'r');
    plot(vmLatsMNeg,(1:length(vmLatsMNeg))./length(vmLatsMNeg),'k');
    t(3) = title('VisMov Cells');
    legend(['Pos: n=',num2str(length(vmLatsMPos))],['Neg: n=',num2str(length(vmLatsMNeg))]);
    set(gca,'YLim',[0 1],'YTick',0:.5:1,'YTickLabel',0:50:100);

    subplot(2,2,4); hold on;
    plot(allPosM,(1:length(allPosM))./length(allPosM),'r');
    plot(allNegM,(1:length(allNegM))./length(allNegM),'k');
    t(4) = title('All Cells');
    legend(['Pos: n=',num2str(length(allPosM))],['Neg: n=',num2str(length(allNegM))]);
    set(gca,'YLim',[0 1],'YTick',0:.5:1,'YTickLabel',0:50:100);
    st      = suptitle(['Pre-Saccade Latencies - Monkey ',upper(monk{im}(1))]);
    sl(1)   = suplabel('Pre-Saccade Latency','x');
    sl(2)   = suplabel('Cumulative Density (%)');
    set(st,'fontsize',16); set(sl,'fontsize',14);
    set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    saveas(gcf,sprintf('%s/%s-PosSpksMovLat.png',saveDir,monk{im}));
    
    
end