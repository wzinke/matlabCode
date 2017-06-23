function klFindExamplev1(varargin)
close all

% Set defaults
area = {'FEF'};
type = 'vis';
monk = 'Gauss';
saccEv = 'SRT';
vRange = -200:800;
mRange = -1000:0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-t','type'}
            type = varargin{varStrInd(iv)+1};
        case {'-a','area'}
            area = varargin{varStrInd(iv)+1};
        case {'-m','monk'}
            monk = varargin{varStrInd(iv)+1};
    end
end

% Set constants
xlFile = './klDataBookKeeping_mg.xlsx';

global excelNum excelAll

% Read in file
[excelNum,~,excelAll] = xlsread(xlFile,monk);

% Get column headers
aCol = find(strcmp(excelAll(4,:),'area'),1);
typeCol = find(strcmp(excelAll(4,:),'typeAlt'),1);

% Find the appropriate rows
myRows = zeros(size(excelAll,1),1);
for ia = 1:length(area),
    myRows = myRows | strcmp(excelAll(:,aCol),area{ia});
end
myRows = find(myRows & strcmp(excelAll(:,typeCol),type));

% Loop through rows
for ir = 1:length(myRows),
    thisRow = myRows(ir);
    [path,file] = klRowToFile(thisRow,'-m',monk);
    load([path,file{1}]);
    
    % Sort spikes by SRT
    mySpks = spiketimes(Task.Correct == 1 & ~isnan(Task.(saccEv)),:);
    mySRT  = Task.(saccEv)(Task.Correct == 1 & ~isnan(Task.(saccEv)));
    if strcmp(saccEv,'SRT'), mySRT = mySRT + Task.GoCue(Task.Correct == 1 & ~isnan(Task.(saccEv))); end
    [sortSacc,sortInd] = sort(mySRT,'ascend');
    mySpks = mySpks(sortInd,:);
    
    vSpks = mySpks;
    mSpks = mySpks-repmat(sortSacc,1,size(mySpks,2));
    for ii = 1:size(sortSacc,1),
        vSpks(ii,vSpks(ii,:) > sortSacc(ii)) = nan;
    end
    [vSDF,vTimes] = klSpkRatev2(vSpks,'-q',1);
    [mSDF,mTimes] = klSpkRatev2(mSpks,'-q',1);
    
    figure(ir);
    %% Plot Vis Response
    subplot(1,2,1);
%     pltMeanStd(vTimes(ismember(vTimes,vRange)),nanmean(vSDF(:,ismember(vTimes,vRange)),1),nanstd(vSDF(:,ismember(vTimes,vRange)),1)./sqrt(sum(~isnan(vSDF(:,ismember(vTimes,vRange))),1)),'k');
    plot(vTimes(ismember(vTimes,vRange)),nanmean(vSDF(:,ismember(vTimes,vRange)),1),'k');
    hold on;
    % Now make the raster
    myRange = get(gca,'YLim'); 
    step = (diff(myRange)./size(vSpks,1))*(2/3);
    for iir = 1:size(vSpks,1),
        plot(vSpks(iir,:),(ones(1,size(vSpks,2)).*myRange(1))+(step*iir),'.','linestyle','none','markersize',1,'color',[.3 .3 .3]);
        plot(sortSacc(iir),myRange(1)+(step*iir),'or','linestyle','none','markersize',1);
    end
    plot(vTimes(ismember(vTimes,vRange)),nanmean(vSDF(:,ismember(vTimes,vRange)),1),'color',[.2 .2 .8],'linewidth',1);
    vl = vline(0); set(vl,'color','k','linewidth',2);
    xlabel('Time From Stim');
    set(gca,'XLim',[vRange(1),vRange(end)]);
    
    %% Now SRT
    subplot(1,2,2);
%     pltMeanStd(mTimes(ismember(mTimes,mRange)),nanmean(mSDF(:,ismember(mTimes,mRange)),1),nanstd(mSDF(:,ismember(mTimes,mRange)),1)./sqrt(sum(~isnan(mSDF(:,ismember(mTimes,mRange))),1)),'k');
    hold on;
    % Now make the raster
    plot(mTimes(ismember(mTimes,mRange)),nanmean(mSDF(:,ismember(mTimes,mRange)),1));
    myRange = get(gca,'YLim'); 
    step = (diff(myRange)./size(mSpks,1))*(2/3);
    for iir = 1:size(mSpks,1),
        plot(mSpks(iir,:),(ones(1,size(mSpks,2)).*myRange(1))+(step*iir),'.','linestyle','none','markersize',1,'color',[.3 .3 .3]);
        plot(-sortSacc(iir),myRange(1)+(step*iir),'og','linestyle','none','markersize',1);
    end
    plot(mTimes(ismember(mTimes,mRange)),nanmean(mSDF(:,ismember(mTimes,mRange)),1),'color',[.2 .2 .8],'linewidth',1);
    vl = vline(0); set(vl,'color','k','linewidth',2);
    set(gca,'XLim',[mRange(1),mRange(end)]);
    xlabel(sprintf('Time From %s',saccEv));
    
    % Get name characteristics
    fName = file{1}(1:11); cellID = file{1}(strfind(file{1},'DSP'):(strfind(file{1},'DSP')+5));
    suptitle(sprintf('%s - %s',fName,cellID));
    
    if mod(ir,10) == 0,
        keyboard;
        close all;
    end
end
keyboard