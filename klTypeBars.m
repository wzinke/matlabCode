%% Set constants
xlFile = 'klDataBookKeeping_mg.xlsx';
monk = {'Gauss','Helmholtz'};

figure();
%% Monkey loop
for im = 1:length(monk)
    %% Load excel file
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    
    %% Get columns
    load xlCols
    vtCol = find(strcmp(excelAll(4,:),'vTrans'));
    vsCol = find(strcmp(excelAll(4,:),'vSust'));
    psCol = find(strcmp(excelAll(4,:),'pSacc'));
    typeCol = find(strcmp(excelAll(4,:),'typeWZ'));
    
    %% Cell types
%     isVis = strcmp(excelAll(:,col.type),'vis');
%     isMov = strcmp(excelAll(:,col.type),'mov');
%     isVM  = strcmp(excelAll(:,col.type),'vismov');
    isVis = strcmp(excelAll(:,typeCol),'vis');
    isMov = strcmp(excelAll(:,typeCol),'mov');
    isVM  = strcmp(excelAll(:,typeCol),'vismov');
    
    %% Get Stat relations
    vtSig = excelNum(:,vtCol) < .05;
    vsSig = excelNum(:,vsCol) < .05;
    psSig = excelNum(:,psCol) < .05 & excelNum(:,psCol+1) >= .05;
    
    %% Get counts and percentages
    % Vis cells
    visVTC  = sum(isVis & vtSig & ~vsSig & ~psSig);
    visVSC  = sum(isVis & ~vtSig & vsSig & ~psSig);
    visVC   = sum(isVis & vtSig & vsSig & ~psSig);
    visPSC  = sum(isVis & ~vtSig & ~vsSig & psSig);
    visVPC   = sum(isVis & (vtSig | vsSig) & psSig);
    
    visVTP  = visVTC/sum(isVis);
    visVSP  = visVSC/sum(isVis);
    visVP   = visVC/sum(isVis);
    visPSP  = visPSC/sum(isVis);
    visVPP  = visVPC/sum(isVis);
    
    % Mov cells
    movVTC  = sum(isMov & vtSig & ~vsSig & ~psSig);
    movVSC  = sum(isMov & ~vtSig & vsSig & ~psSig);
    movVC   = sum(isMov & vtSig & vsSig & ~psSig);
    movPSC  = sum(isMov & ~vtSig & ~vsSig & psSig);
    movVPC   = sum(isMov & (vtSig | vsSig) & psSig);
    
    movVTP  = movVTC/sum(isMov);
    movVSP  = movVSC/sum(isMov);
    movVP   = movVC/sum(isMov);
    movPSP  = movPSC/sum(isMov);
    movVPP  = movVPC/sum(isMov);
    
    % VM cells
    vmVTC  = sum(isVM & vtSig & ~vsSig & ~psSig);
    vmVSC  = sum(isVM & ~vtSig & vsSig & ~psSig);
    vmVC   = sum(isVM & vtSig & vsSig & ~psSig);
    vmPSC  = sum(isVM & ~vtSig & ~vsSig & psSig);
    vmVPC   = sum(isVM & (vtSig | vsSig) & psSig);
    
    vmVTP  = vmVTC/sum(isVM);
    vmVSP  = vmVSC/sum(isVM);
    vmVP   = vmVC/sum(isVM);
    vmPSP  = vmPSC/sum(isVM);
    vmVPP  = vmVPC/sum(isVM);
    
    % Make bars
    barData = [visVTP, movVTP, vmVTP;...
               visVSP, movVSP, vmVSP;...
               visVP,  movVP,  vmVP;...
               visPSP, movPSP, vmPSP;...
               visVPP, movVPP, vmVPP];
    
%     figure();
    subplot(1,2,im); hold on; % Stacked
    bar(barData','stacked');
    legend('VT','VS','VT&S','PS','V&P');
    
%     figure(2); subplot(1,2,im); hold on;
%     bar(barData');
%     legend('VT','VS','VT&S','PS','V&P');
%     
end
   
suptitle(excelAll(4,typeCol));    
    