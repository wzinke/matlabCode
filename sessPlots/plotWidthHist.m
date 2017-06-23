%% Start timing the script

clearvars
close all;

startTic = tic;

fprintf('\n\n******************************************\n');
fprintf('******     STARTING SHELL ANALYSIS    *****\n');
fprintf('******************************************\n\n\n');

%% Set user defined variables (when/if this becomes a function, these will be set in "varargin")
monk        = {'Gauss','Helmholtz'};                     % Cell array of strings for monkey names to analyze. Even if it's one monkey, still keep it a cell
type        = 'DSP';                             % DSP for spikes, LFP for LFP, and EEG for EEG
task        = 'Capture';
physType    = 'waves';
upsamp      = 10;                               % For spline interpolations, by what factor should time values be upsampled
outFold     = 'C:\Users\Kaleb\Documents\MATLAB\Plots\widthHists';

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
widthBins = 0:32;
runAreas = {};

%% Set and move to root matlab folder for this computer
% Set
matlabFold = 'C:/Users/Kaleb/Documents/MATLAB';
% Move
if ~strcmp(cd,matlabFold), cd(matlabFold); end

% Make Wolf's scripts available
addpath(genpath('./wz_code'));

%% Let's get the data in the mix
dataFold = 'Y:/Users/Wolf/ephys_db';

%% Start monkey loop
fprintf('**** BEGINNING MONKEY LOOP ****\n\n');
for im = 1:length(monk)
    fprintf('Fetching data from %s...\n',monk{im});
    clear sessFolders numSess goodSess axH
    %close all;
    sessFolders = dir(sprintf('%s/%s/20*',dataFold,monk{im}));
    numSess = length(sessFolders);
    
    gridLocs = cell(1,numSess);
    goodSess = 1;
    wvWidth = cell(24,1,2);
    
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
            % the M-L and A-P coordinates and area
            if sum(isnan(unitNum)) == length(unitNum)
                gridLocs{is} = [ml_coor,ap_coor];
                areas{is}        = area(~ismember(area,'?'));
                if ~ismember(areas{is},runAreas),
                    runAreas{length(runAreas)+1} = areas{is};
                end
            end
            arInd = find(ismember(runAreas,areas{is}));
            if arInd > size(wvWidth,2)
                wvWidth = cat(2,wvWidth,cell(size(wvWidth,1),(arInd - size(wvWidth,2)),2));
            end
            % Get info on this unit
            unitNum(iu)     = abc2num(wave.suffix);
            chanNum(iu)     = str2double(DSPname(4:end-1));
            if isnan(unitNum(iu)) || isnan(chanNum(iu)),
                keyboard
            end
            maxWV(iu)       = max(max(wave.waves));
            [wvMean(iu,:), wvStd(iu,:)] = getMeanWaves(wave.waves,'-pm',0,'median');
            
            % Check for positive spikes
            isPos(iu)       = abs(max(wvMean(iu,:))) > abs(min(wvMean(iu,:)));
            isPosThresh(iu) = wave.thresh > 0;
%             if isPosThresh(iu) && ~isPos(iu)
%                 keyboard
%             end
            
            %% Get waveform width
            wvWidth{chanNum(iu),arInd,isPos(iu)+1}(length(wvWidth{chanNum(iu),arInd,isPos(iu)+1})+1) = wv_width(wvMean(iu,:),'-t',wvTimes);
                        
        end
        %% End unit loop (iu)
       
        %keyboard
    end
    %% End session loop (is)
     
    %% Make and plot histograms of wave width
    %x = figure();
    % Make the histograms
    allHists = cell(size(wvWidth));
    wbCell   = cell(size(wvWidth)); [wbCell{:}] = deal(widthBins);
    allHists = cellfun(@hist,wvWidth,wbCell,'Uniformoutput',0);
    %maxHist = max(max(max(cellfun(@max,allHists))));
    maxHist = max(max(cellfun(@max,allHists),[],3),[],1);
    for ia = 1:length(runAreas)
        figure(ia);
        axH = depthAxes('-fwx',ia,.35,.45*(im-1));
        for ip = 1:2
            for ic = 1:length(axH)
                if ic == 1, t = title(axH(ic),sprintf('Monkey %s',upper(monk{im}(1)))); set(t,'fontsize',16); end
                barH  = bar(axH(ic),allHists{ic,ia,ip}); set(barH,'facecolor',colors(ip)); 
                ylabel(axH(ic),num2str(ic));
                set(axH(ic),'YLim',[0 max(maxHist(ia),.5)],'XTick',widthBins(1):4:widthBins(end));
            end
        end
        a(1) = xlabel('Spike Width');
        set(axH(end),'XTickLabel',widthBins(1):4:widthBins(end));
        if im == 1, [ax, a(2)] = suplabel('Channel #','y'); end; clear ax
        set(a,'fontsize',14,'fontweight','bold');
        if im == length(monk)
            st = suptitle(sprintf('Area %s',runAreas{ia})); set(st,'fontsize',20,'fontweight','bold');
            set(gcf,'paperposition',[.2 .1 10.5 7.5],'papersize',[11 8]);
            areaTitle = runAreas{ia}; areaTitle(ismember(areaTitle,'/\')) = '-';
            saveas(gcf,sprintf('%s/%s-WidthHistograms.bmp',outFold,areaTitle));
        end
        
    end
    
end
%% End monkey loop (im)
fprintf('***  Script completed in %s   ***\n',printTiming(startTic));




