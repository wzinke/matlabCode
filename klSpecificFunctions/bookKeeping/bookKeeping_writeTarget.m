clearvars;
close all;

%% Start timing the script
startTic = tic;

fprintf('\n\n******************************************\n');
fprintf('******     STARTING SHELL ANALYSIS    *****\n');
fprintf('******************************************\n\n\n');

%% Set user defined variables (when/if this becomes a function, these will be set in "varargin")
monk        = {'Helmholtz','Gauss'};                     % Cell array of strings for monkey names to analyze. Even if it's one monkey, still keep it a cell
type        = 'DSP';                             % DSP for spikes, LFP for LFP, and EEG for EEG
task        = 'Capture';
physType    = 'rates';

% Print out user defined variables to command line
fprintf('---> Analyzing %s: %s data from task %s, monkey %s',physType,type,task,monk{1});
for im = 2:length(monk)
    fprintf(', %s',monk{im});
end
fprintf('\n\n');

%% Initialize some variables that will be used later
thisSub = 0;
outFile = 'cosDataBookKeeping_capture.xlsx';

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
            rf = getRF(spiketimes,Task);
            [pVals, pNames] = klGetTarget(spiketimes,Task);
            
            outRow = {thisSess,sessUnits(iu).name,chanNum,area,'',isPosThresh,isPos,''};
            suffLoop = {'Targ','Loc 0','Loc 90','Loc 180','Loc 270'};
            startColTxt = {'I','N','S','X','AC'};
            startCol = cellfun(@abc2num,startColTxt);
            for il = 1:length(suffLoop),
                stInd = find(~cellfun(@isempty,strfind(pNames,suffLoop{il})),1);
                for ip = 0:3
                    outRow = cat(2,outRow,pVals(stInd+ip));
                end
                outRow{length(outRow)+1} = '';
            end
            outMat = cat(1,outMat,outRow);
            
        end
    end
    xlswrite(outFile,outMat,monk{im},'A4');
%     keyboard
end
