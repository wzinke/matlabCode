%% Start timing the script

clearvars

startTic = tic;

fprintf('\n\n******************************************\n');
fprintf('******     STARTING SHELL ANALYSIS    *****\n');
fprintf('******************************************\n\n\n');

%% Set user defined variables (when/if this becomes a function, these will be set in "varargin")
%monk        = {'Gauss','Helmholtz'};                     % Cell array of strings for monkey names to analyze. Even if it's one monkey, still keep it a cell
monk = {'Helmholtz'};
type        = 'DSP';                             % DSP for spikes, LFP for LFP, and EEG for EEG
task        = 'MG';
physType    = 'waves';
upsamp      = 10;                               % For spline interpolations, by what factor should time values be upsampled

funType     = 'pltMeanWaves';

% Print out user defined variables to command line
fprintf('---> Analyzing %s: %s data from task %s, monkey %s',physType,type,task,monk{1});
for im = 2:length(monk)
    fprintf(', %s',monk{im});
end
fprintf('\n\n');

%% Initialize some variables that will be used later
thisSub = 0;
colors = 'rgbcmk';
runAreas = {};
areaTitles = {};

%% Set and move to root matlab folder for this computer
% Set
matlabFold = 'C:/Users/Kaleb/Documents/MATLAB';
% Move
if ~strcmp(cd,matlabFold), cd(matlabFold); end

% Make scripts from other folders available
addpath(genpath('./wz_code'));
addpath(genpath('./externFunctions'));
addpath(genpath('./klFunctions'));

%% Let's get the data in the mix
dataFold = 'Y:/Users/Wolf/ephys_db';

%% Start monkey loop
fprintf('**** BEGINNING MONKEY LOOP ****\n\n');
for im = 1:length(monk)
    fprintf('Fetching data from %s...\n',monk{im});
    clear sessFolders numSess goodSess a
    close all;
    sessFolders = dir(sprintf('%s/%s/20*',dataFold,monk{im}));
    numSess = length(sessFolders);
    whichFig = [];
    whichSub = [];

    gridLocs = cell(1,numSess);
    goodSess = 1;
    %% Start session loop
    for is = 1:numSess
        clear thisSess sessUnits wvMean wvStd
        fprintf('\tAnalyzing Session %s (%d/%d)...',sessFolders(is).name,is,numSess);
        
        thisSess    = sessFolders(is).name;
        sessUnits   = dir(sprintf('%s/%s/%s/%s/%s*',dataFold,monk{im},thisSess,type,type));
        chanNum     = nan(length(sessUnits),1);
        unitNum     = nan(length(sessUnits),1);
        
        fprintf('Found %d units\n',length(sessUnits));
        
        wvMean      = nan(length(sessUnits),32);
        wvStd       = nan(length(sessUnits),32);
        maxWV       = nan(length(sessUnits),1);
        wvTimes     = 1:32;
        
        %% Start unit loop
        for iu = 1:length(sessUnits)
            clear f fName
            % Check for valid data files
            f = dir(sprintf('%s/%s/%s/%s/%s/%s_%s_%s*',dataFold,monk{im},thisSess,type,sessUnits(iu).name,thisSess,sessUnits(iu).name,task));
            if isempty(f),
                fprintf('\t\tNo "%s" file found for unit %s... Skipping to next unit\n',task,sessUnits(iu).name);
                continue;
            elseif length(f) > 1,
                fprintf('\t\tMultiple "%s" files found for unit %s... Please check directory structure\n',task,sessUnits(iu).name);
                fprintf('\t\t\tDefaulting to first file in list\n');
                %keyboard
            end
            fName = f(1).name;
                
            % Load unit data
            loadTic = tic;
            switch physType
                case 'rates'
                    load(sprintf('%s/%s/%s/%s/%s/%s',dataFold,monk{im},thisSess,type,sessUnits(iu).name,fName));
                case 'waves'
                    load(sprintf('%s/%s/%s/%s/%s/waves/%s_waves.mat',dataFold,monk{im},thisSess,type,sessUnits(iu).name,fName(1:(end-4))));
            end
            % Next line can be used to time this load process, but clutters
            % the command line
            %fprintf('\t\tUnit data loaded in %s\n',printTiming(loadTic));
            
            %% This section was to allow me to play around with some functions to learn what they do
            %  I'm commenting this out for now, but it'll likely return
%             switch funType
%                 case 'SPK_resp_profile'
%                     figure(iu);
%                     RESP = SPK_resp_profile(spiketimes,Task.SaccEnd,0,1);
%                 case 'pltMeanWaves'
%                     [wvMean, wvStd] = pltMeanWaves(wave.waves,'-ctp','b',wave.thresh,0);
%                     if wave.thresh > 0
%                         fprintf('Found a positive-spiking cell!\n');
%                         keyboard
%                     end
%                     %case
%             end
            %% End learning section
            
            %% Get waveform info
            
            % If this is the first unit in this session/penetration, get
            % the M-L and A-P coordinates
            if sum(isnan(unitNum)) == length(unitNum)
                gridLocs{is} = [ml_coor,ap_coor];
                areas{is}        = area(~ismember(area,'?'));
                if ~ismember(areas{is},runAreas),
                    runAreas{length(runAreas)+1} = areas{is};
                    thisTitle = areas{is}; thisTitle(ismember(thisTitle,'/\')) = '-';
                    areaTitles{length(areaTitles)+1} = thisTitle;
                end
            end
            arInd = find(ismember(runAreas,areas{is}));
            if arInd > size(whichFig,2)
                whichFig = cat(2,whichFig,nan(1,arInd - size(whichFig,2)));
                whichSub = cat(2,whichSub,nan(1,arInd - size(whichSub,2)));
            end
            
            unitNum(iu)     = abc2num(wave.suffix);
            chanNum(iu)     = str2double(DSPname(4:end-1));
            if isnan(unitNum(iu)) || isnan(chanNum(iu)),
                keyboard
            end
            maxWV(iu)       = max(max(wave.waves));
            thresh(iu)      = wave.thresh;
            [wvMean(iu,:), wvStd(iu,:)] = getMeanWaves(wave.waves,'-pm',0,'median');
            
            %% Get waveform width
            wvWidth(iu) = wv_width(wvMean(iu,:),'-t',wvTimes);
                        
        end
        %% End unit loop (iu)
        
        %% Plot average waveforms from units, when present
        % Open new figure, if needed (5 days per fig)
        if sum(isnan(maxWV(:))) == numel(maxWV),
            continue
        end
        yDelt = ceil(max(maxWV)/2);
        
        if isnan(whichFig(arInd)), whichFig(arInd) = 1; end
        if isnan(whichSub(arInd)), whichSub(arInd) = 1; end;
        thisSub = whichSub(arInd); figOff = 1;
        while thisSub > 5, thisSub = thisSub - 5; figOff = figOff + 1; end
        figure(figOff+(10*(arInd-1))); subplot(1,5,thisSub); hold on;
        %if mod(goodSess,5) == 1, thisFig = figure(); thisSub = 1; else thisSub = thisSub + 1; end
        %subplot(1,5,thisSub); hold on;
        for iu = 1:length(sessUnits)
            if ~isnan(chanNum(iu)) && ~isnan(unitNum(iu))
                pltMeanStd(wvTimes,wvMean(iu,:)-(yDelt*chanNum(iu)),wvStd(iu,:),colors(unitNum(iu)));
                h(1) = hline(-(yDelt*chanNum(iu))); set(h(1),'linewidth',1,'color','k');
                h(2) = hline(thresh(iu)-(yDelt*chanNum(iu))); set(h(2),'linewidth',1,'color','k','linestyle','--');
                clear h;
            end
        end
        set(gca,'YLim',[-(max(chanNum)*yDelt)-ceil(max(maxWV)),0+ceil(max(maxWV))],'YTick',(-yDelt*max(chanNum)):yDelt:(-yDelt),'YTickLabel',max(chanNum):-1:1);
        t = title(sprintf('%s - %s',area,thisSess)); set(t,'fontsize',16,'fontweight','bold');
        a(1) = xlabel('Time Bin');
        if thisSub == 1, a(2) = ylabel('Channel Number'); end; 
        set(a,'fontsize',14,'fontweight','bold');
        goodSess = goodSess + 1;
        whichSub(arInd) = whichSub(arInd)+1;
        %keyboard
    end
    %% End session loop (is)
    
    % Save figures created in this subject loop
    openFigs = cell2mat(get(get(0,'Children'),'Number'));
    set(openFigs,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
    for iFig = 1:length(openFigs)
        saveas(openFigs(iFig),sprintf('./Plots/Waveforms/%s-%s-WFplots-SessionSet-%d.bmp',monk{im},areaTitles{ceil(openFigs(iFig)/10)},mod(openFigs(iFig),10)));
    end
end
%% End monkey loop (im)
fprintf('***  Script completed in %s   ***\n',printTiming(startTic));




