function [outMean, outStd, outTimes, outEvents] = klPlotChansv3(excelRow,varargin)

% Start the clock
startTic = tic;

%% Add script folders to path
if ~any(ismember(path,'./wz_code'));
    addpath(genpath('./wz_code'));
end
addpath(genpath('./klFunctions'));
addpath(genpath('./externFunctions'));

%% Set defaults
alignEvents     = {'StimOnset','SRT'};
events          = {'StimOnset','SRT','RewardTone','Reward'};
%colors          = 'rgbcmkrgbcmkrgbcmkrgbcmkrgbcmkrgbcmkrgbcmk';
colors          = repmat([.8 .2 .2; .2 .8 .2; .2 .2 .8;.2 .8 .8; .8 .2 .8; 0 0 0],5,1);
xlFile          = './klDataBookKeeping_mg.xlsx';
monk            = 'Helmholtz';
avWvType        = {'Mean','Median'};
offset          = 0;
dataFold        = 'Y:/Users/Wolf/ephys_db';
task            = 'MG';
doWZ            = 1;
doSave          = 0;
evntWind        = 100;
baseWind        = 100;
visSustWind     = [150, 250];
saveDir         = './Plots/UnitSummary';
quiet           = 0;
doZ             = 0;
clip            = 0;
clipEv          = 'SRT';
doLat           = 1;

%% Decode varargin
% Check that there are enough input arguments, if optional flags are
% set
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case '-a'
            alignEvents = varargin{varStrInd(iv)+1};
        case '-e'    
            events      = varargin{varStrInd(iv)+1};
        case '-m'
            monk        = varargin{varStrInd(iv)+1};
        case '-o'
            offset      = varargin{varStrInd(iv)+1};
        case '-f'
            figNum      = varargin{varStrInd(iv)+1};
%             case 'w'
%                 doWZ        = varargin{iv};
        case '-s'
            doSave      = varargin{varStrInd(iv)+1};
        case '-w'
            evntWind    = varargin{varStrInd(iv)+1};
        case '-b'
            baseWind    = varargin{varStrInd(iv)+1};
        case '-q'
            quiet       = varargin{varStrInd(iv)+1};
        case '-z'
            doZ         = varargin{varStrInd(iv)+1};
        case '-c'
            clip        = varargin{varStrInd(iv)+1};
    end
end

%% Define axes positions based on aligned events (max 2 supported now)
switch length(alignEvents)
    case 1
        rastAx = {[.10 .10 .80 .40]};
        rateAx = {[.10 .60 .80 .15]};
        waveAx = {[.10 .85 .80 .15]};
    case 2
        rastAx = {[.10 .10 .35 .40], [.55 .10 .35 .40]};
        rateAx = {[.10 .60 .35 .15], [.55 .60 .35 .15]};
        waveAx = {[.10 .85 .35 .15], [.55 .85 .35 .15]};
       
end

%% Load in the excel file data and make dimensions consistent
[excelNum,excelText,excelAll] = xlsread(xlFile,monk);
% excelNum = cat(1,nan(1,size(excelNum,2)),excelNum);
% excelNum = cat(2,nan(size(excelNum,1),size(excelAll,2)-size(excelNum,2)),excelNum);

hRow    = find(strcmp(excelAll(:,2),'Name'));

outMean = cell(length(excelRow),2);
outTimes = cell(length(excelRow),2);
outStd  = cell(length(excelRow),2);
outEvents = cell(length(excelRow),2);

for ir = 1:length(excelRow)
    clear evMat alMat
    fprintf('Processing Channel %d (of %d - %d%%)\n',ir,length(excelRow),floor(ir*100/length(excelRow)));
    %% Determine which row to plot and get basic information from excel file
    thisRow     = excelRow(ir)-offset;
    thisSess    = excelAll{thisRow,strcmp(excelAll(hRow,:),'Name')};
    thisChan    = excelAll{thisRow,strcmp(excelAll(hRow,:),'chanCode')};
    if doWZ
        pVis        = excelNum(thisRow,strcmp(excelAll(hRow,:),'pVis(WZ)'));
        pMov        = excelNum(thisRow,strcmp(excelAll(hRow,:),'pMov(WZ)'));
    else
        pVis        = excelNum(thisRow,strcmp(excelAll(hRow,:),'vTrans'));
        pMov        = excelNum(thisRow,strcmp(excelAll(hRow,:),'pSacc'));
    end
    
    % Load requested data file
    f = dir(sprintf('%s/%s/%s/DSP/%s/%s_%s_%s*',dataFold,monk,thisSess,thisChan,thisSess,thisChan,task));
    if length(f) == 0,
        if ~quiet
            fprintf('\tNo "%s" file found for unit %s... Skipping to next unit\n',task,thisChan);
        end
        continue;
    elseif length(f) > 1,
        if ~quiet,
            fprintf('\tMultiple "%s" files found for unit %s... Please check directory structure\n',task,thisChan);
        end
        %keyboard
    end
    fName = f(1).name(1:end-4);
    
    load(sprintf('%s/%s/%s/DSP/%s/%s',dataFold,monk,thisSess,thisChan,fName));
    load(sprintf('%s/%s/%s/DSP/%s/waves/%s_waves',dataFold,monk,thisSess,thisChan,fName));
    
    if exist('area','var'), areaTitle = eval('area'); end
    areaTitle(ismember(areaTitle,'/\')) = '-';
    areaTitle(ismember(areaTitle,'?'))  = '';
    
	for ie = 1:length(alignEvents),
		alMat(:,ie) = [Task.(alignEvents{ie})];
	end
    for ie = 1:length(events)
        evMat(:,ie) = [Task.(events{ie})];
    end
    corrTrials = Task.Correct;
    
%     if clip,
%         spiketimes(spiketimes > repmat(Task.(clipEv),1,size(spiketimes,2))) = nan;
%     end
%     
    for ia = 1:length(alignEvents)
        clear blWind tstWind blSpks tstSpks
        %% Get the alignment times and adjust the spike times accordingly
        alignTimes = Task.(alignEvents{ia});
        alignSpikes = spiketimes - repmat(alignTimes,1,size(spiketimes,2));
        
        %alMat = alMat - repmat(alignTimes,1,size(alMat,2));
        evMat = evMat - repmat(alignTimes,1,size(evMat,2));
        
% 		alignTimes = alignTimes(valTrials);
%         alignSpikes = alignSpikes(valTrials);
%         evMat = evMat(valTrials,:);
% 		
        
        % Get SDF/time window
%         [sessSDF, tRange] = klSpkRatev2(alignSpikes,'-q',1);
        % Plot data in appropriate axes
        if ~exist('figNum','var'); figure(); figNum = get(gcf,'Number'); end
        figure(figNum+(ir-1));
        rast(ia) = axes('position',rastAx{ia});
        rate(ia) = axes('position',rateAx{ia});
        wvax(ia) = axes('position',waveAx{ia});
        
        % Sort the trials according to the other (of the two) events
        if ia == 1,
            [sortVals, sortInds] = sort(alMat(:,2)-alignTimes,1,'ascend');
        elseif ia == 2
            [sortVals, sortInds] = sort(alMat(:,1)-alignTimes,1,'descend');
        end
        %sortSDF = sessSDF(sortInds,:);
        sortSpikes = alignSpikes(sortInds,:);
        sortAlMat = alMat(sortInds,:);
        sortEvMat = evMat(sortInds,:);
        sortTimes = alignTimes(sortInds,:);
        sortCorr  = corrTrials(sortInds);
        
        % Limit to just valid trials (valid = non-nan for both alignment
        % times
        %valTrials = ~isnan(sortAlMat(:,1)) & ~isnan(sortAlMat(:,2));
        valTrials = sortCorr == 1;
        %valSDF    = sortSDF(valTrials,:);
        valSpikes = sortSpikes(valTrials,:);
        valTimes  = sortTimes(valTrials);
        valEvMat  = sortEvMat(valTrials,:);
        
        [valSDF,tRange] = klSpkRatev2(valSpikes,'-q',1);
        tVals           = tRange;
        
        if clip
            if ~ismember(alignEvents{ia},{'SRT','SaccEnd'}),
                valSpikes(valSpikes > repmat(valEvMat(:,2),1,size(valSpikes,2))) = nan;
                for it = 1:size(valEvMat,1),
                    lastTimeInd = find(tVals > valEvMat(it,2),1);
    %                 if it == 63,
    %                     keyboard
    %                 end
                    if ~isempty(lastTimeInd)
                        valSDF(it,lastTimeInd:end) = nan;
                    end
               end
                tVals = tVals(~isnan(nanmean(valSDF,1)));
            else
                valSpikes(valSpikes < repmat(valEvMat(:,1)-baseWind,1,size(valSpikes,2))) = nan;
                for it = 1:size(valEvMat,1),
                    firstTimeInd = find(tVals < (valEvMat(it,1)-baseWind),1,'last');
    %                 if it == 63,
    %                     keyboard
    %                 end
                    if ~isempty(firstTimeInd)
                        valSDF(it,1:firstTimeInd) = nan;
                    end
               end
                %tVals = tVals(~isnan(nanmean(valSDF,1)));
            end
        end
        
        
        plotRange = [min([-150,floor(min(valEvMat(:,1))/100)*100]), ceil(max(Task.Reward-alignTimes)/100)*100];
        tInds = (tVals >= plotRange(1) & tVals <= plotRange(2));
        
        %% Take out nans in valSDF (Still need to figure out how these nans get there...)
%         tVals(isnan(nanmean(valSDF,1))) = [];
%         valSDF(:,isnan(nanmean(valSDF,1))) = [];
%         
        
        %% Plot Rasters
        axes(rast(ia)); hold on;
        [blSpks, tstSpks] = deal(nan(1,size(valSpikes,1)));
        [blWind, tstWind, vsWind] = deal(nan(size(valSpikes,1),2));
        for it = 1:size(valSpikes,1)
            % Goal: turn this whole section into klPlotRast and put the
            % below looped sections (calculating the patches) into another
            % it loop...
            
            % Plot this trial raster
            hp = plot(valSpikes(it,:),ones(1,size(valSpikes,2)).*it,'.k','linestyle','none','markersize',3,'markerfacecolor','k');
            
            % Get this trial event windows
            blWind(it,:)        = [valEvMat(it,1) - baseWind, valEvMat(it,1)];
            tstWind(it,:)       = [valEvMat(it,ia) - (evntWind*mod(ia+1,2)), valEvMat(it,ia) + (evntWind*mod(ia,2))];
            vsWind(it,:)        = [valEvMat(it,1) + visSustWind(1), valEvMat(it,1) + visSustWind(2)];
            
            % Get spikes in those windows
            blSpks(it)          = (sum(valSpikes(it,:) >= blWind(it,1)  & valSpikes(it,:) < blWind(it,2)))/diff(blWind(it,:)./1000);
            tstSpks(it)         = (sum(valSpikes(it,:) >= tstWind(it,1) & valSpikes(it,:) < tstWind(it,2)))/diff(tstWind(it,:)./1000);
            
        end
        pVal = signrank(blSpks,tstSpks);
        for ie = 1:length(events)
            %plot(sortEvMat(valTrials,ie)-valTimes,1:length(valTimes),'o','color',colors(ie,:),'markersize',1,'markerfacecolor',colors(ie,:));
            plot(sortEvMat(valTrials,ie),1:length(valTimes),'o','color',colors(ie,:),'markersize',1,'markerfacecolor',colors(ie,:));
        end
        blPatch  = patch([blWind(:,1); flipud(blWind(:,2))], [1:sum(valTrials), sum(valTrials):-1:1],'c');
        tstPatch = patch([tstWind(:,1);flipud(tstWind(:,2))],[1:sum(valTrials),sum(valTrials):-1:1],'m');
        vsPatch  = patch([vsWind(:,1);flipud(vsWind(:,2))],[1:sum(valTrials),sum(valTrials):-1:1],'g');
        
        set([blPatch,tstPatch,vsPatch],'EdgeAlpha',0,'FaceAlpha',.3);
        
        set(gca,'XLim',plotRange,'YLim',[1, length(valTimes)]);
        a(1) = xlabel(sprintf('Time (aligned on %s)',alignEvents{ia}));
        if ia == 1, a(2) = ylabel('Trial Number'); end
        set(a,'fontsize',12);
        %keyboard
        
        blMean = nanmean(blSpks); blStd = nanstd(blSpks);
        
        %% Add values to output variables
        if doZ, outData = normData(valSDF,'-z',[blMean, blStd]); else outData = valSDF; end
        outMean{ir,ia}      = nanmean(outData,1);
        outTimes{ir,ia}     = (floor(tRange(1)):ceil(tRange(2)));
        outStd{ir,ia}       = nanstd(outData,1)./sqrt(size(outData,1));
        outEvents{ir,ia}    = nanmean(valEvMat,1);
        
        %% Plot SDF
        axes(rate(ia));
        %pltMeanStd((floor(tRange(1)):ceil(tRange(2))),nanmean(valSDF,1),nanstd(valSDF,1)./sqrt(size(valSDF,1)),'k');
        pltMeanStd(tVals(tInds),nanmean(valSDF(:,tInds),1),nanstd(valSDF(:,tInds),1)./sqrt(sum(isfinite(valSDF(:,tInds)),1)),'k');
        for ie = 1:length(events)
            vl(ia,ie) = vline(nanmean(sortEvMat(valTrials,ie))); set(vl(ia,ie),'color',colors(ie,:),'linestyle','-','linewidth',3);
        end
        if ia == 1, ylabel('Firing Rate (Hz)','fontsize',12); end
        t = title(sprintf('p = %.3f',pVal)); set(t,'fontsize',14);
        yLim = get(gca,'YLim');
        
        blRatePatch(ia) = patch([nanmean(blWind(:,1)),nanmean(blWind(:,1)),nanmean(blWind(:,2)),nanmean(blWind(:,2))],[yLim,fliplr(yLim)],'c');
        set(blRatePatch(ia),'EdgeAlpha',0,'FaceAlpha',.3);
        tstRatePatch(ia) = patch([nanmean(tstWind(:,1)),nanmean(tstWind(:,1)),nanmean(tstWind(:,2)),nanmean(tstWind(:,2))],[yLim,fliplr(yLim)],'m');
        set(tstRatePatch(ia),'EdgeAlpha',0,'FaceAlpha',.3);
        vsRatePatch(ia)  = patch([nanmean(vsWind(:,1)),nanmean(vsWind(:,1)),nanmean(vsWind(:,2)),nanmean(vsWind(:,2))],[yLim,fliplr(yLim)],'g');
        set(vsRatePatch(ia),'EdgeAlpha',0,'FaceAlpha',.3);
        set(gca,'XLim',plotRange);
        
        if doLat,
            if ia == 1,
                [thisLat, latTimes] = klPoissLatv2(alignSpikes,'-rwd',Task.Reward-alignTimes,'-stim',Task.StimOnset-alignTimes);
            elseif ia == 2,
                [thisLat, mLatTimes] = klPoissLatv2(alignSpikes,'-rwd',Task.Reward-alignTimes,'-stim',Task.StimOnset-alignTimes,'-stop',Task.SRT-alignTimes,'direction','rev');
            end
            latLine = vline(thisLat); set(latLine,'color','k','linestyle','--','linewidth',2);
%             keyboard
        end
        
        
        %% Plot Waveform
        axes(wvax(ia));
        [wvMean, ~] = getMeanWaves(wave.waves,'-pm',1,avWvType{ia});
        hZ = hline(0); set(hZ,'linewidth',3);
        hT = hline(wave.thresh); set(hT,'linewidth',2,'linestyle','--');
        if ia == 1, ylabel('WF (uV)'); end;
        t = title(sprintf('%s WF',avWvType{ia})); 
        
        if abs(max(wvMean)) > abs(min(wvMean)),
            wfSign = 'Positive';
        elseif abs(max(wvMean)) < abs(min(wvMean)),
            wfSign = 'Negative';
        else
            if wave.thresh > 0
                wfSign = 'Positive';
            else
                wfSign = 'Negative';
            end
        end
        
        
    end
    % Check which rate axis has the greater yLim
    if length(alignEvents) > 1,
        yLims = cell2mat(get(rate,'YLim'));
        set(rate,'YLim',yLims(find(yLims(:,2) == max(yLims(:,2)),1),:));
        set(tstRatePatch,'YData',[yLims(find(yLims(:,2) == max(yLims(:,2)),1),:),fliplr(yLims(find(yLims(:,2) == max(yLims(:,2)),1),:))]);
        set(blRatePatch,'YData',[yLims(find(yLims(:,2) == max(yLims(:,2)),1),:),fliplr(yLims(find(yLims(:,2) == max(yLims(:,2)),1),:))]);
        set(vsRatePatch,'YData',[yLims(find(yLims(:,2) == max(yLims(:,2)),1),:),fliplr(yLims(find(yLims(:,2) == max(yLims(:,2)),1),:))]);
        set(vl,'YData',yLims(find(yLims(:,2) == max(yLims(:,2)),1),:));
    end
    
    st = suptitle(sprintf('%s - %s (%s) - %s',monk,thisSess,eval('area'),thisChan)); set(st,'fontsize',16);
    %keyboard
    if doSave
        set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(gcf,sprintf('%s/%s/%s-%s-%s-%s.png',saveDir,wfSign,monk,areaTitle,thisSess,thisChan));
        close(gcf);
    end
end
%keyboard
fprintf('***  Script completed in %s   ***\n',printTiming(startTic));
% winopen(saveDir);
