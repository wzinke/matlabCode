% klPlotSessions

%% Set changeable defaults
fresh = 0;

areas = {'FEF'};
task = 'MG';
doRF = 0;
zDim = 2;
zType = 'baseline';
vTimes = -200:300;
mTimes = -300:200;
wvTMin = -200; wvTMax = 500;

%% Set constants info
monk = {'Gauss','Helmholtz','Darwin'};
xlFile = 'klDataBookKeeping_mg.xlsx';
root = 'Y:/Users/Wolf/ephys_db';
getCNum = @(chanCode) str2double(chanCode(4:5));
saveDir = './sessPlots';

%% Set globals
global excelNum excelAll

%% Start monkey loop
for im = 1:length(monk),
    % Load in xlFile
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    
    %% Get sessions from the right area
    areaCol = find(strcmp(excelAll(4,:),'area'));
    uSess = unique(excelAll(find(ismember(excelAll(5:end,areaCol),areas))+4,2));
    
    % Loop through sessions
    for is = 1:length(uSess),
        fprintf('Doing monkey %s: Session %d of %d\n',im,is,length(uSess));
        close all;
        thisSess = uSess{is};
        if ~fresh && exist(sprintf('%s/%s-%s-Norm.pdf',saveDir,monk{im},thisSess),'file'),
            continue
        end
        
        %% Get units from this session
        sessRows = find(strcmp(excelAll(:,2),thisSess));
%         sessChans = cellfun(@getCNum,excelAll(sessRows,strcmp(excelAll(4,:),'chanCode')));
        
        [sessSDF, sessTimes, sessLFP, sessLFPTimes, sessWaves, sessWvTimes, normMatUnit, normMatLFP] = klGetFilev5(sessRows,'-m',monk{im},'-t',task,'-rf',doRF,'-z',zDim,'ztype',zType);
        
        sessChans = nan(length(sessRows),1);
        for iu = 1:length(sessRows),
            sessChans(iu) = getCNum(excelAll{sessRows(iu),strcmp(excelAll(4,:),'chanCode')});
        end
        
%         visAx=depthAxes('-w',.15,'-f',is,'m',.05);
%         visLFPAx = depthAxes('-w',.15,'x',.2,'-f',is,'m',.05);
%         movAx=depthAxes('-w',.15,'x',.4,'-f',gcf,'m',.05);
%         movLFPAx = depthAxes('-w',.15,'x',.6,'-f',is,'m',.05);
%         wvAx = depthAxes('-w',.1,'x',.8,'-f',gcf,'m',.05);
        clear visMap visLFPMap movMap movLFPMap
        for iu = 1:length(sessRows),
            visMap(iu,1:length(vTimes))     = sessSDF{iu,1}(ismember(sessTimes{iu,1},vTimes));
            visLFPMap(iu,1:length(vTimes))  = sessLFP{iu,1}(ismember(sessLFPTimes{iu,1},vTimes));
            movMap(iu,1:length(mTimes))     = sessSDF{iu,2}(ismember(sessTimes{iu,2},mTimes));
            movLFPMap(iu,1:length(mTimes))  = sessLFP{iu,2}(ismember(sessLFPTimes{iu,2},mTimes));
        end
        
        visMapNorm = (visMap-repmat(min(visMap,[],2),1,size(visMap,2)))./(repmat(max(visMap,[],2)-min(visMap,[],2),1,size(visMap,2)));
        visLFPNorm = (visLFPMap-repmat(min(visLFPMap,[],2),1,size(visLFPMap,2)))./(repmat(max(visLFPMap,[],2)-min(visLFPMap,[],2),1,size(visLFPMap,2)));
        movMapNorm = (movMap-repmat(min(movMap,[],2),1,size(movMap,2)))./(repmat(max(movMap,[],2)-min(movMap,[],2),1,size(movMap,2)));
        movLFPNorm = (movLFPMap-repmat(min(movLFPMap,[],2),1,size(movLFPMap,2)))./(repmat(max(movLFPMap,[],2)-min(movLFPMap,[],2),1,size(movLFPMap,2)));
        
        visMapZ = (visMap-repmat(normMatUnit(:,1),1,size(visMap,2)))./(repmat(normMatUnit(:,2),1,size(visMap,2)));
        visLFPZ = (visLFPMap-repmat(normMatLFP(:,1),1,size(visLFPMap,2)))./(repmat(normMatLFP(:,2),1,size(visLFPMap,2)));
        movMapZ = (movMap-repmat(normMatUnit(:,1),1,size(movMap,2)))./(repmat(normMatUnit(:,2),1,size(movMap,2)));
        movLFPZ = (movLFPMap-repmat(normMatLFP(:,1),1,size(movLFPMap,2)))./(repmat(normMatLFP(:,2),1,size(movLFPMap,2)));
        
        figure(is);
        subplot(1,4,1);
        imagesc(visMap); colorbar;
        set(gca,'XTick',find(vTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('Unit-Vis'); ylabel('Unit Count');
        
        subplot(1,4,2);
        imagesc(visLFPMap); colorbar;
        set(gca,'XTick',find(vTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('LFP-Vis');
        
        subplot(1,4,3);
        imagesc(movMap); colorbar;
        set(gca,'XTick',find(mTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('Unit-Vis');
        
        subplot(1,4,4);
        imagesc(movLFPMap); colorbar;
        set(gca,'XTick',find(mTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('LFP-Sacc');
        
        suptitle(sprintf('%s - %s - Raw',monk{im},thisSess));
        set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(gcf,sprintf('%s/%s-%s-Raw.pdf',saveDir,monk{im},thisSess));
        
        figure(is+100);
        subplot(1,4,1);
        imagesc(visMapNorm); colorbar;
        set(gca,'XTick',find(vTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('Unit-Vis'); ylabel('Unit Count');
        
        subplot(1,4,2);
        imagesc(visLFPNorm); colorbar;
        set(gca,'XTick',find(vTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('LFP-Vis');
        
        subplot(1,4,3);
        imagesc(movMapNorm); colorbar;
        set(gca,'XTick',find(mTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('Unit-Vis');
        
        subplot(1,4,4);
        imagesc(movLFPNorm); colorbar;
        set(gca,'XTick',find(mTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('LFP-Sacc');
        
        suptitle(sprintf('%s - %s - Norm',monk{im},thisSess));
        set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(gcf,sprintf('%s/%s-%s-Norm.pdf',saveDir,monk{im},thisSess));
        
        figure(is+200);
        subplot(1,4,1);
        imagesc(visMapZ); colorbar;
        set(gca,'XTick',find(vTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('Unit-Vis'); ylabel('Unit Count');
        
        subplot(1,4,2);
        imagesc(visLFPZ); colorbar;
        set(gca,'XTick',find(vTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('LFP-Vis');
        
        subplot(1,4,3);
        imagesc(movMapZ); colorbar;
        set(gca,'XTick',find(mTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('Unit-Vis');
        
        subplot(1,4,4);
        imagesc(movLFPZ); colorbar;
        set(gca,'XTick',find(mTimes==0),'XTickLabel',0,'tickdir','out');
        xlabel('LFP-Sacc');
        
        suptitle(sprintf('%s - %s - Z-Scored',monk{im},thisSess));
        set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(gcf,sprintf('%s/%s-%s-Z.pdf',saveDir,monk{im},thisSess));
        
    end
end
        
        
        
        

