clearvars;
close all;

%% Start timing the script
startTic = tic;

fprintf('\n\n******************************************\n');
fprintf('******     STARTING SHELL ANALYSIS    *****\n');
fprintf('******************************************\n\n\n');

%% Set user defined variables (when/if this becomes a function, these will be set in "varargin")
monk        = {'Gauss','Helmholtz'};                     % Cell array of strings for monkey names to analyze. Even if it's one monkey, still keep it a cell
type        = 'DSP';                             % DSP for spikes, LFP for LFP, and EEG for EEG
task        = 'MG';
physType    = 'rates';

funType     = 'pltMeanWaves';

% Print out user defined variables to command line
fprintf('---> Analyzing %s: %s data from task %s, monkey %s',physType,type,task,monk{1});
for im = 2:length(monk)
    fprintf(', %s',monk{im});
end
fprintf('\n\n');

%% Initialize some variables that will be used later
thisSub = 0;
outFile = 'klDataBookKeeping_mg.xlsx';

%% Set and move to root matlab folder for this computer
% Set
matlabFold = 'C:/Users/Kaleb/Documents/MATLAB';
% Move
if ~strcmp(cd,matlabFold), cd(matlabFold); end

% Make Wolf's scripts available
addpath(genpath('./wz_code'));

%% Let's get the data in the mix
dataFold = 'Y:/Users/Wolf/ephys_db';

chanLoaded = 0;
unitTypes = {};

%% Start monkey loop
fprintf('**** BEGINNING MONKEY LOOP ****\n\n');
for im = 1:length(monk)
    fprintf('Fetching data from %s...\n',monk{im});
    clear sessFolders numSess outMat
    sessFolders = dir(sprintf('%s/%s/20*',dataFold,monk{im}));
    numSess = length(sessFolders);
    outMat = {};
    %% Start session loop
    for is = 1:numSess
        clear thisSess sessUnits
        fprintf('\tAnalyzing Session %s...',sessFolders(is).name);
        thisSess = sessFolders(is).name;
        sessUnits = dir(sprintf('%s/%s/%s/%s/%s*',dataFold,monk{im},thisSess,type,type));
        fprintf('Found %d units\n',length(sessUnits));
        
        %% Start unit loop
        for iu = 1:length(sessUnits)
            clear f fName
            % Check for valid data files
            f = dir(sprintf('%s/%s/%s/%s/%s/%s_%s_%s*',dataFold,monk{im},thisSess,type,sessUnits(iu).name,thisSess,sessUnits(iu).name,task));
            if length(f) == 0,
                fprintf('\tNo "%s" file found for unit %s... Skipping to next unit\n',task,sessUnits(iu).name);
                continue;
            elseif length(f) > 1,
                fprintf('\tMultiple "%s" files found for unit %s... Please check directory structure\n',task,sessUnits(iu).name);
                %keyboard
            end
            fName = f(1).name;
                
            % Load unit data
            loadTic = tic;
            
            load(sprintf('%s/%s/%s/%s/%s/%s',dataFold,monk{im},thisSess,type,sessUnits(iu).name,fName));
            load(sprintf('%s/%s/%s/%s/%s/waves/%s_waves.mat',dataFold,monk{im},thisSess,type,sessUnits(iu).name,fName(1:(end-4))));
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
            chanLoaded = chanLoaded + 1;
            
            % Check if the area string needs changed
            if sum(isletter(area)) == 0
                area = [' ',area];
            end
            area(strfind(area,'?')) = [];
            
            chanNum     = str2double(DSPname(4:end-1));
            [wvMean, wvStd] = getMeanWaves(wave.waves,'-pm',0,'median');
            
            % Check for positive spikes
            isPos       = abs(max(wvMean)) > abs(min(wvMean));
            isPosThresh = wave.thresh > 0;
         
            % Get waveform width
            wvWidth = wv_width(wvMean,'-t',1:length(wvMean));
            
            % Get response property classification
            [pVals, pTypes, unitTypes{length(unitTypes)+1}, vmi, numSpks] = klGetType(spiketimes,Task);
            resp = SPK_resp_profile(spiketimes,Task.SaccEnd);
            
            % Get latency
            %lat = poissLat_Long(spiketimes,Task);
            [sdf, sdfTime] = klSpkRatev2(spiketimes);
            lat = klGetLat(sdf,'-b',[-100, 0]+nanmedian(Task.StimOnset),'-f',5,'-t',sdfTime,'-type','ttest','-filt',Task.Correct == 1);
            
            
            
            %% Write to excel
            %         Session   Channel             Depth    Area  Type    
            %         VM Index  pVis                         pMov 
            %         Type(WZ)       pVis(WZ)         pMov(WZ)         Positive Spike(WF)?
            %         Positive Spike (T)?  spkWidth
            
           
            outRow = {thisSess ,sessUnits(iu).name, chanNum, area, unitTypes{length(unitTypes)}, resp.resptype, vmi, '',...
                pVals(strcmp(pTypes,'vTrans')), pVals(strcmp(pTypes,'vSust')), pVals(strcmp(pTypes,'pSacc')), pVals(strcmp(pTypes,'VSvsPreSacc')), pVals(strcmp(pTypes,'postSacc')),...
                pVals(strcmp(pTypes,'VTvsPostSacc')), pVals(strcmp(pTypes,'rwdTone')), pVals(strcmp(pTypes,'Rwd')), pVals(strcmp(pTypes,'preErr')), pVals(strcmp(pTypes,'Err')),...
                '', double(isPosThresh), double(isPos), wvWidth};
            outMat = cat(1,outMat,outRow); 
            %keyboard
        end
        %keyboard
        %% End unit loop (iu)
    end
    %% End session loop (is)
    
    xlswrite(outFile,outMat,monk{im},'A5');
end
%% End monkey loop (im)
fprintf('***  Script completed in %s   ***\n',printTiming(startTic));


winopen(outFile);

