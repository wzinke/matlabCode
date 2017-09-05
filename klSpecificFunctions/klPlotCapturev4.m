function [outMean, outStd, outTimes, outEvents] = klPlotCapturev4(excelRow,varargin)

% Start the clock
startTic = tic;

%% Add script folders to path
if ~any(ismember(path,'./wz_code'));
    addpath(genpath('./wz_code'));
end
addpath(genpath('./klFunctions'));
addpath(genpath('./externFunctions'));

%% Set Constants
colors          = repmat([.8 .2 .2; .2 .8 .2; .2 .2 .8;.2 .8 .8; .8 .2 .8; 0 0 0],5,1);
xlFile          = './klDataBookKeeping_mg.xlsx';
avWvType        = {'Mean','Median'};
dataFold        = 'Y:/Users/Wolf/ephys_db';
task            = 'Capture';
tInCol          = 'r';
tOutCol         = 'k';
degreeToSubplotFour = [6, 2, 4, 8];
degreeToSubplotEight = [6 3 2 1 4 7 8 9];
saveDir         = './Plots/unitSummary/targetSelection';
inRF            = 180;
outRF           = 0;
nStdsY           = 7;

%% Set defaults
alignEvents     = {'StimOnset','SRT'};
events          = {'StimOnset','SRT','RewardTone','Reward'};
monk            = 'Helmholtz';
offset          = 0;
doWZ            = 0;
doSave          = 0;
evntWind        = 100;
baseWind        = 100;
visSustWind     = [150, 300];
quiet           = 0;
doZ             = 0;
clip            = 1;
clipEv          = 'SRT';
doStats         = 1;

%% Decode varargin
if nargin > 1
    % Check that there are enough input arguments, if optional flags are
    % set
    if length(varargin{1}) ~= length(varargin), error('Not enough input arguments'); end
    for iv = 2:length(varargin{1})
        switch varargin{1}(iv)
            case 'a'
                alignEvents = varargin{iv};
            case 'e'    
                events      = varargin{iv};
            case 'm'
                monk        = varargin{iv};
            case 'o'
                offset      = varargin{iv};
            case 'f'
                figNum      = varargin{iv};
%             case 'w'
%                 doWZ        = varargin{iv};
            case 's'
                doSave      = varargin{iv};
            case 'w'
                evntWind    = varargin{iv};
            case 'b'
                baseWind    = varargin{iv};
            case 'v'
                visSustWind = varargin{iv};
            case 'q'
                quiet       = varargin{iv};
            case 'z'
                doZ         = varargin{iv};
            case 'c'
                clip        = varargin{iv};
            case 't'
                doStats     = varargin{iv};
        end
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
% excelNum = cat(1,nan(size(excelAll,1)-size(excelNum,1),size(excelNum,2)),excelNum);
% excelNum = cat(2,nan(size(excelNum,1),size(excelAll,2)-size(excelNum,2)),excelNum);

hRow = find(strcmp(excelAll(:,2),'Name'));

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
        pVis        = excelNum(thisRow,strcmp(excelAll(hRow,:),'pVis'));
        pMov        = excelNum(thisRow,strcmp(excelAll(hRow,:),'pMov'));
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
    maxLocMn  = 0; maxLocStd = 0; minLocMn = inf; minLocStd = inf;
    for ia = 1:length(alignEvents)
        clear blWind tstWind blSpks tstSpks
        maxTypeMn = 0; maxTypeStd = 0; minTypeMn = inf; minTypeStd = inf;
        %% Get the alignment times and adjust the spike times accordingly
        alignTimes = Task.(alignEvents{ia});
        alignSpikes = spiketimes - repmat(alignTimes,1,size(spiketimes,2));
        
        evMat = evMat - repmat(alignTimes,1,size(evMat,2));
        
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
        sortSRT   = Task.SRT(sortInds);
        sortTLoc  = Task.TargetLoc(sortInds);
        if isfield(Task,'DistLoc'), sortDLoc = Task.DistLoc(sortInds); else sortDLoc = nan(1,length(sortInds)); end;
        if isfield(Task,'Singleton'), sortSing = Task.Singleton(sortInds); else sortSing = nan(1,length(sortInds)); end;
        
        % Limit to just valid trials (valid = non-nan for both alignment
        % times
        %valTrials = ~isnan(sortAlMat(:,1)) & ~isnan(sortAlMat(:,2));
        valAl     = ~isnan(sortAlMat(:,1));
        for iai = 2:size(sortAlMat,2),
            valAl = valAl & ~isnan(sortAlMat(:,iai));
        end
        valTrials = sortCorr == 1 & valAl;
%         valTrials = sortCorr == 1;
        valSpikes = sortSpikes(valTrials,:);
        valTimes  = sortTimes(valTrials);
        valEvMat  = sortEvMat(valTrials,:);
        valSRT    = sortSRT(valTrials);
        valTLoc   = sortTLoc(valTrials,:);
        valDLoc = sortDLoc(valTrials);
        valSing = sortSing(valTrials);
        
        [valSDF,tRange] = klSpkRatev2(valSpikes,'-q',1);
        tVals           = floor(tRange(1)):ceil(tRange(end));
        
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
        
        valSTD    = nanstd(valSDF,1)./sqrt(sum(isfinite(valSDF),1));
        maxSTD    = max(valSTD(sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75)));
        maxMn     = max(valSDF(:));
        plotRange = [min([-150,floor(min(sortEvMat)/100)*100]), ceil(max(Task.RewardTone-alignTimes)/100)*100];
        tInds = (tVals >= plotRange(1) & tVals <= plotRange(2));
        
            
        %% Take out nans in valSDF (Still need to figure out how these nans get there...)
%         tVals(isnan(nanmean(valSDF,1))) = [];
%         valSDF(:,isnan(nanmean(valSDF,1))) = [];
%         
        tInInd      = find(valTLoc == inRF);
        tOutInd     = find(valTLoc == outRF);
        allTLocs    = unique(valTLoc(isfinite(valTLoc)));
        if length(allTLocs) == 4,
            degreeToSubplot = degreeToSubplotFour;
            %continue
        elseif length(allTLocs) == 8
            degreeToSubplot = degreeToSubplotEight;
            %keyboard
        end
        
        dSalInInd   = find(valDLoc == inRF);
        dSalOutInd  = find(valDLoc == outRF);
        
        %% Plot Rasters
        axes(rast(ia)); hold on;
        [blSpks, tstSpks, vsSpks] = deal(nan(1,size(valSpikes,1)));
        [blWind, tstWind, vsWind] = deal(nan(size(valSpikes,1),2));
        [aovTLoc, aovDLoc, aovSing]  = deal(nan(1,size(valSpikes,1)));
        numTPlotted = 0;
        [aovTargLoc, aovDistLoc, aovSLoc, aovTrialSpikes] = deal([]);
        for il = 1:length(allTLocs)
            thisTInd = find(valTLoc == allTLocs(il));
            aovTLoc(thisTInd) = repmat(il,1,length(thisTInd));
            for it = 1:length(thisTInd)
                % Plot this trial raster
                numTPlotted = numTPlotted + 1;
                hp = plot(valSpikes(thisTInd(it),:),ones(1,size(valSpikes,2)).*numTPlotted,'.','color',colors(il,:),'linestyle','none','markersize',3,'markerfacecolor',colors(il,:));

                % Get this trial event windows
                blWind(thisTInd(it),:)        = [valEvMat(thisTInd(it),1) - baseWind, valEvMat(thisTInd(it),1)];
                tstWind(thisTInd(it),:)       = [valEvMat(thisTInd(it),ia) - (evntWind*mod(ia+1,2)), valEvMat(thisTInd(it),ia) + (evntWind*mod(ia,2))];
                vsWind(thisTInd(it),:)        = [valEvMat(thisTInd(it),1) + visSustWind(1), valEvMat(thisTInd(it),1) + visSustWind(2)];

                % Get spikes in those windows
                blSpks(thisTInd(it))          = (sum(valSpikes(thisTInd(it),:) >= blWind(thisTInd(it),1)  & valSpikes(thisTInd(it),:) < blWind(thisTInd(it),2)))/diff(blWind(thisTInd(it),:)./1000);
                tstSpks(thisTInd(it))         = (sum(valSpikes(thisTInd(it),:) >= tstWind(thisTInd(it),1) & valSpikes(thisTInd(it),:) < tstWind(thisTInd(it),2)))/diff(tstWind(thisTInd(it),:)./1000);
                vsSpks(thisTInd(it))          = (sum(valSpikes(thisTInd(it),:) >= vsWind(thisTInd(it),1) & valSpikes(thisTInd(it),:) < vsWind(thisTInd(it),2)))/diff(vsWind(thisTInd(it),:)./1000);
                aovDLoc(thisTInd(it))         = valDLoc(thisTInd(it));
                aovSing(thisTInd(it))         = valSing(thisTInd(it)); 
                
                aovTargLoc     = cat(2,aovTargLoc,valTLoc(thisTInd(it)));
                aovDistLoc     = cat(2,aovDistLoc,valDLoc(thisTInd(it)));
                aovSLoc        = cat(2,aovSLoc, valSing(thisTInd(it)));
                aovTrialSpikes = cat(2,aovTrialSpikes, tstSpks(thisTInd(it)));
                
            end
        end
        
        pVal = anovan(tstSpks,{aovTLoc},'display','off');
        numTPlotted = 0;
        clear blPatch tstPatch vsPatch
        for il = 1:length(allTLocs)
            thisTInd = find(valTLoc == allTLocs(il));
            for ie = 1:length(events)
                %plot(sortEvMat(thisTInd,ie),1+numTPlotted:length(thisTInd)+numTPlotted,'o','color',colors(ie,:),'markersize',1,'markerfacecolor',colors(ie,:));
                plot(valEvMat(thisTInd,ie),1+numTPlotted:length(thisTInd)+numTPlotted,'o','color',colors(ie,:),'markersize',1,'markerfacecolor',colors(ie,:));
            end
            blPatch(il)  = patch([blWind(thisTInd,1); flipud(blWind(thisTInd,2))], [1+numTPlotted:length(thisTInd)+numTPlotted, length(thisTInd)+numTPlotted:-1:1+numTPlotted],'c');
            tstPatch(il) = patch([tstWind(thisTInd,1);flipud(tstWind(thisTInd,2))],[1+numTPlotted:length(thisTInd)+numTPlotted,length(thisTInd)+numTPlotted:-1:1+numTPlotted],'m');
            %vsPatch(il)  = patch([vsWind(thisTInd,1);flipud(vsWind(thisTInd,2))],[1+numTPlotted:length(thisTInd)+numTPlotted,length(thisTInd)+numTPlotted:-1:1+numTPlotted],'g');
            
            numTPlotted = numTPlotted + length(thisTInd);
        end
        set(rast(ia),'YLim',[0 size(valSDF,1)]);
        %set([blPatch,tstPatch,vsPatch],'EdgeAlpha',0,'FaceAlpha',.3);
        set([blPatch,tstPatch],'EdgeAlpha',0,'FaceAlpha',.3);
        
        clear a
        %set(gca,'XLim',[floor(tRange(1)),ceil(tRange(2))],'YLim',[1, length(valTimes)]);
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
        axes(rate(ia)); hold on;
        %pltMeanStd((floor(tRange(1)):ceil(tRange(2))),nanmean(valSDF,1),nanstd(valSDF,1)./sqrt(size(valSDF,1)),'k');
        clear checkRF
        for il = 1:length(allTLocs),
            thisTInd = find(valTLoc == allTLocs(il));
            pltMeanStd(tVals(tInds),nanmean(valSDF(thisTInd,tInds),1),nanstd(valSDF(thisTInd,tInds),[],1)./sqrt(sum(isfinite(valSDF(thisTInd,tInds)),1)),colors(il,:));
            checkRF(il) = nanmean(nanmean(valSDF(thisTInd,tInds),1),2);
            
            theseCols  = sum(isfinite(valSDF(thisTInd,:)),1) >= (max(sum(isfinite(valSDF(thisTInd,:)),1))*.75);
%             maxLocMn  = max([maxLocMn,nanmean(valSDF(thisTInd,sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75)),1)]);
%             maxLocStd = max([maxLocStd,nanstd(valSDF(thisTInd,sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75)),1)./sqrt(sum(isfinite(valSDF(:,sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75))),1))]);
%             minLocMn  = min([minLocMn,nanmean(valSDF(thisTInd,sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75)),1)]);
%             minLocStd = min([minLocStd,nanstd(valSDF(thisTInd,sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75)),1)./sqrt(sum(isfinite(valSDF(:,sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75))),1))]);
            maxLocMn  = max([maxLocMn,nanmean(valSDF(thisTInd,theseCols),1)]);
            maxLocStd = max([maxLocStd,nanstd(valSDF(thisTInd,theseCols),1)./sqrt(sum(isfinite(valSDF(thisTInd,theseCols)),1))]);
            minLocMn  = min([minLocMn,nanmean(valSDF(thisTInd,theseCols),1)]);
            minLocStd = min([minLocStd,nanstd(valSDF(thisTInd,theseCols),1)./sqrt(sum(isfinite(valSDF(thisTInd,theseCols)),1))]);
            
        end
        thisRF(ir,ia)       = allTLocs(checkRF == max(checkRF));
        thisAntiRF(ir,ia)   = allTLocs(checkRF == min(checkRF));
        thisLatTarg     = excelNum(thisRow,find(strcmp(excelAll(hRow,:),'vLat-Target'))+(ia-1));
        thisLatType     = excelNum(thisRow,find(strcmp(excelAll(hRow,:),'vLat-Type'))+(ia-1));
    
        for ie = 1:length(events)
            vl(ia,ie) = vline(nanmedian(sortEvMat(valTrials,ie))); set(vl(ia,ie),'color',colors(ie,:),'linestyle','-','linewidth',3);
            
        end
        if doStats
            [~, sigTimes] = klGetLatencyv4(spiketimes, Task, Task.TargetLoc, '-rav',thisRF(ir,ia),{alignEvents{ia}},0);
            if ismember(alignEvents{ia},{'SRT','SaccEnd'}),
                pScat     = scatter(sigTimes(sigTimes > nanmedian(valEvMat(:,1))),ones(1,sum(sigTimes > nanmedian(valEvMat(:,1)))).*(maxLocMn+(nStdsY-2)*maxLocStd),'*r');
            else
                pScat     = scatter(sigTimes(sigTimes < nanmedian(valEvMat(:,2))),ones(1,sum(sigTimes < nanmedian(valEvMat(:,2)))).*(maxLocMn+(nStdsY-2)*maxLocStd),'*r');
            end
        end
        latL(ia) = vline(thisLatTarg); set(latL(ia),'color','k','linestyle',':','linewidth',3);
        
        %set(gca,'XLim',[floor(tRange(1)),ceil(tRange(2))]);
        if ia == 1, ylabel('Firing Rate (Hz)','fontsize',12); end
        t = title(sprintf('p = %.3f',pVal)); set(t,'fontsize',14);
        yLim = get(gca,'YLim');
        
        set([rate(ia),rast(ia)],'XLim',plotRange);
        
        blRatePatch(ia) = patch([nanmedian(blWind(:,1)),nanmedian(blWind(:,1)),nanmedian(blWind(:,2)),nanmedian(blWind(:,2))],[yLim,fliplr(yLim)],'c');
        set(blRatePatch(ia),'EdgeAlpha',0,'FaceAlpha',.3);
        tstRatePatch(ia) = patch([nanmedian(tstWind(:,1)),nanmedian(tstWind(:,1)),nanmedian(tstWind(:,2)),nanmedian(tstWind(:,2))],[yLim,fliplr(yLim)],'m');
        set(tstRatePatch(ia),'EdgeAlpha',0,'FaceAlpha',.3);
%         vsRatePatch(ia)  = patch([nanmedian(vsWind(:,1)),nanmedian(vsWind(:,1)),nanmedian(vsWind(:,2)),nanmedian(vsWind(:,2))],[yLim,fliplr(yLim)],'g');
%         set(vsRatePatch(ia),'EdgeAlpha',0,'FaceAlpha',.3);
        
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
        
        
        %% Now that we have the RF, we can plot Target In, Salient Distractor (Singleton), and Non-Salient Distractor in RF
        %  Get parts of valSDF where T is in  the different locations, Sal.
        %  Distractor in that location, Nonsalient distractor
        figure(((figNum+(ir-1))*100) + ia); 
        for il = 1:length(allTLocs)
            clear aovAllSpks thisTargSpks thisSalSpks thisNSalSpks
            subplot(3,3,degreeToSubplot(il)); hold on;
            valPoints    = tInds; %numTrialsSDF > 0;
            
            
            pltMeanStd(tVals(valPoints),nanmean(valSDF(valTLoc == allTLocs(il),valPoints),1),nanstd(valSDF(valTLoc == allTLocs(il),valPoints),1)./sqrt(sum(isfinite(valSDF(valTLoc == allTLocs(il),valPoints)),1)),[.8, .2, .2]);
            pltMeanStd(tVals(valPoints),nanmean(valSDF(valDLoc == allTLocs(il) & valSing == 1,valPoints),1),nanstd(valSDF(valDLoc == allTLocs(il) & valSing == 1,valPoints),1)./sqrt(sum(isfinite(valSDF(valDLoc == allTLocs(il) & valSing == 1,valPoints)),1)),[.2, .8, .2]);
            pltMeanStd(tVals(valPoints),nanmean(valSDF(valDLoc == allTLocs(il) & valSing == 0,valPoints),1),nanstd(valSDF(valDLoc == allTLocs(il) & valSing == 0,valPoints),1)./sqrt(sum(isfinite(valSDF(valDLoc == allTLocs(il) & valSing == 0,valPoints)),1)),[.2 .2 .8]);
            
            trialCrit = {valTLoc == allTLocs(il),valDLoc == allTLocs(il) & valSing == 1,valDLoc == allTLocs(il) & valSing == 0};
            for it = 1:3
%                 maxTypeMn  = max([maxTypeMn,nanmean(valSDF(trialCrit{it},sum(isfinite(valSDF(trialCrit{it},:)),1) >= (size(valSDF(trialCrit{it},:),1)*.75)),1)]);
%                 maxTypeStd = max([maxTypeStd,nanstd(valSDF(trialCrit{it},sum(isfinite(valSDF(trialCrit{it},:)),1) >= (size(valSDF(trialCrit{it},:),1)*.75)),1)./sqrt(sum(isfinite(valSDF(:,sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75))),1))]);
%                 minTypeMn  = min([minTypeMn,nanmean(valSDF(trialCrit{it},sum(isfinite(valSDF(trialCrit{it},:)),1) >= (size(valSDF(trialCrit{it},:),1)*.75)),1)]);
%                 minTypeStd = min([minTypeStd,nanstd(valSDF(trialCrit{it},sum(isfinite(valSDF(trialCrit{it},:)),1) >= (size(valSDF(trialCrit{it},:),1)*.75)),1)./sqrt(sum(isfinite(valSDF(:,sum(isfinite(valSDF),1) >= (size(valSDF,1)*.75))),1))]);
               % dimCrit    = sum(isfinite(valSDF(trialCrit{it},:
                theseCols  = sum(isfinite(valSDF(trialCrit{it},:)),1) >= (max(sum(isfinite(valSDF(trialCrit{it},:)),1))*.75);
                maxTypeMn  = max([maxTypeMn,nanmean(valSDF(trialCrit{it},theseCols),1)]);
                maxTypeStd = max([maxTypeStd,nanstd(valSDF(trialCrit{it},theseCols),1)./sqrt(sum(isfinite(valSDF(trialCrit{it},theseCols)),1))]);
                minTypeMn  = min([minTypeMn,nanmean(valSDF(trialCrit{it},theseCols),1)]);
                minTypeStd = max([minTypeStd,nanstd(valSDF(trialCrit{it},theseCols),1)./sqrt(sum(isfinite(valSDF(trialCrit{it},theseCols)),1))]);
                
            end
            
            thisTargSpks = aovTrialSpikes(aovTargLoc == allTLocs(il));
            thisSalSpks = aovTrialSpikes(aovDistLoc == allTLocs(il) & aovSLoc == 1);
            thisNSalSpks = aovTrialSpikes(aovDistLoc == allTLocs(il) & aovSLoc == 0);
            
            aovAllSpks = cat(2,thisTargSpks,thisSalSpks,thisNSalSpks);
            aovAllType = cat(2,ones(1,length(thisTargSpks)),ones(1,length(thisSalSpks)).*2,ones(1,length(thisNSalSpks)).*3);
            [pVal(il), t, stats] = anovan(aovAllSpks,{aovAllType},'display','off');
            
            yVals(il,:) = get(gca,'YLim');
            if clip
                if ~ismember(alignEvents{ia},{'SRT','SaccEnd'}),
                    set(gca,'XLim',[plotRange(1),nanmedian(Task.SRT - alignTimes)]);
                else
                    set(gca,'XLim',[nanmedian(blWind(:,1)),plotRange(2)]);
                end
            else
                set(gca,'XLim',plotRange);
            end
            %keyboard
            
        end
        pValAdj = pAdjust(pVal);
        clear blTargPatch tstTargPatch
        clear vlTarg
        for il = 1:length(allTLocs),
            subplot(3,3,degreeToSubplot(il));
%             set(gca,'YLim',[min(yVals(:,1)),max(yVals(:,2))]);
%             set(gca,'YLim',[max([0,minLocMn-10*maxLocStd]),(maxTypeMn+10*maxTypeStd)]); 
            set(gca,'YLim',[max([0,minTypeMn-nStdsY*maxTypeStd]),(maxTypeMn+nStdsY*maxTypeStd)]); 
            %blTargPatch(il) = patch([nanmedian(blWind(:,1)),nanmedian(blWind(:,1)),nanmedian(blWind(:,2)),nanmedian(blWind(:,2))],[min(yVals(:,1)),max(yVals(:,2)),max(yVals(:,2)),min(yVals(:,1))],'c');
            %tstTargPatch(il) = patch([nanmedian(tstWind(:,1)),nanmedian(tstWind(:,1)),nanmedian(tstWind(:,2)),nanmedian(tstWind(:,2))],[min(yVals(:,1)),max(yVals(:,2)),max(yVals(:,2)),min(yVals(:,1))],'m');
            blTargPatch(il) = patch([nanmedian(blWind(:,1)),nanmedian(blWind(:,1)),nanmedian(blWind(:,2)),nanmedian(blWind(:,2))],[max([0,minLocMn-10*maxLocStd]),(maxTypeMn+10*maxTypeStd),fliplr([max([0,minLocMn-10*maxLocStd]),(maxTypeMn+10*maxTypeStd)])],'c');
            tstTargPatch(il) = patch([nanmedian(tstWind(:,1)),nanmedian(tstWind(:,1)),nanmedian(tstWind(:,2)),nanmedian(tstWind(:,2))],[max([0,minLocMn-10*maxLocStd]),(maxTypeMn+10*maxTypeStd),fliplr([max([0,minLocMn-10*maxLocStd]),(maxTypeMn+10*maxTypeStd)])],'m');
            vlTarg(il) = vline(0);
            if pValAdj(il) < .05
                latLType(il) = vline(thisLatType); set(latLType(il),'color','k','linestyle',':','linewidth',3);
            end
            if allTLocs(il) == thisRF(ir,ia),
                t = title(sprintf('RF - p = %.3f',pValAdj(il)));
            elseif allTLocs(il) == thisAntiRF(ir,ia),
                t = title(sprintf('Anti-RF - p = %.3f',pValAdj(il)));
            else
                t = title(sprintf('Neither - p = %.3f',pValAdj(il)));
            end
            
            %keyboard
        end
        set([blTargPatch,tstTargPatch],'facealpha',.3,'edgealpha',0);
        set(vlTarg,'color',colors(1,:));
        
        subplot(3,3,5); hold on;
        [wvMean, ~] = getMeanWaves(wave.waves,'-pm',1,avWvType{ia});
        hZ = hline(0); set(hZ,'linewidth',3);
        hT = hline(wave.thresh); set(hT,'linewidth',2,'linestyle','--');
        if ia == 1, ylabel('WF (uV)'); end;
        t = title(sprintf('%s WF',avWvType{ia})); 
        
        st = suptitle(sprintf('%s - %s (%s) - %s',monk,thisSess,eval('area'),thisChan)); set(st,'fontsize',16);
        a(1) = suplabel('Firing Rate (Hz)','y');
        a(2) = suplabel(sprintf('Time (aligned on %s)',alignEvents{ia}),'x');
        set(a,'fontsize',14,'fontweight','bold');
        if doSave, 
            set(((figNum+(ir-1))*100) + ia,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
            saveas(((figNum+(ir-1))*100) + ia,sprintf('%s/%s/%s/%s-%s-%s-%s.png',saveDir,wfSign,alignEvents{ia},monk,areaTitle,thisSess,thisChan));
            close(((figNum+(ir-1))*100) + ia);
        end

        %keyboard
    end
    
    figure(figNum+(ir-1));
    % Check which rate axis has the greater yLim
    if length(alignEvents) > 1,
        yLims = cell2mat(get(rate,'YLim'));
%         set(rate,'YLim',yLims(find(yLims(:,2) == max(yLims(:,2)),1),:));
        set(rate,'YLim',[max([0,minLocMn-nStdsY*maxLocStd]),maxLocMn+nStdsY*maxLocStd]);
%         set(tstRatePatch,'YData',[yLims(find(yLims(:,2) == max(yLims(:,2)),1),:),fliplr(yLims(find(yLims(:,2) == max(yLims(:,2)),1),:))]);
%         set(blRatePatch,'YData',[yLims(find(yLims(:,2) == max(yLims(:,2)),1),:),fliplr(yLims(find(yLims(:,2) == max(yLims(:,2)),1),:))]);
        %set(vsRatePatch,'YData',[yLims(find(yLims(:,2) == max(yLims(:,2)),1),:),fliplr(yLims(find(yLims(:,2) == max(yLims(:,2)),1),:))]);
%         set(vl,'YData',yLims(find(yLims(:,2) == max(yLims(:,2)),1),:));
        set([tstRatePatch,blRatePatch],'YData',[minLocMn-nStdsY*maxLocStd,maxLocMn+nStdsY*maxLocStd,maxLocMn+nStdsY*maxLocStd,minLocMn-nStdsY*maxLocStd]);
        set(vl,'YData',[max([0,minLocMn-nStdsY*maxLocStd]),maxLocMn+nStdsY*maxLocStd]);
        if exist('latL','var'),
            set(latL,'YData',[max([0,minLocMn-nStdsY*maxLocStd]),maxLocMn+nStdsY*maxLocStd]);
        end
    end
    
    st = suptitle(sprintf('%s - %s (%s) - %s',monk,thisSess,eval('area'),thisChan)); set(st,'fontsize',16);
    %keyboard
    if doSave
        set(figNum+(ir-1),'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
        saveas(figNum+(ir-1),sprintf('%s/%s/%s-%s-%s-%s.png',saveDir,wfSign,monk,areaTitle,thisSess,thisChan));
        close(figNum+(ir-1));
    end
end
%keyboard
fprintf('***  Script completed in %s   ***\n',printTiming(startTic));
% winopen(saveDir);
