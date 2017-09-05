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
task        = 'Capture';
physType    = 'rates';
statSampRate = 5;
visualize = 0;

% Set cutoffs for unreasonable TST values
tstMax = 300;

% Print out user defined variables to command line
fprintf('---> Analyzing %s: %s data from task %s, monkey %s',physType,type,task,monk{1});
for im = 2:length(monk)
    fprintf(', %s',monk{im});
end
fprintf('\n\n');

%% Initialize some variables that will be used later
thisSub = 0;
mgFile  = 'klDataBookKeeping_mg.xlsx';
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
    outMat = {};
    oldSess = '';
    % Load in the excel file data and make dimensions consistent
    [excelNum,excelText,excelAll] = xlsread(mgFile,monk{im});
    writeCol = find(strcmp(excelAll(4,:),'vLat-Target'),1);
    
%     fprintf('Row Loop starting at %s\n',printTiming(startTic));
        
    for ir = 5:size(excelNum,1),
        thisRow     = ir;
        thisSess    = excelAll{thisRow,strcmp(excelAll(4,:),'Name')};
        thisChan    = excelAll{thisRow,strcmp(excelAll(4,:),'chanCode')};
        thisType    = excelAll{thisRow,strcmp(excelAll(4,:),'typeAlt')};
        
        outRow = cell(1,30); [outRow{:}] = deal('');
        
        if ~strcmp(thisSess,oldSess),
            fprintf('\tAnalyzing Session %s...',thisSess);
            chanThisSess = sum(strcmp(excelAll(:,1),thisSess));
            fprintf('Found %d channels\n',chanThisSess);
            oldSess = thisSess;
            chanSoFar = 0;
        end
        chanSoFar = chanSoFar + 1;
        if chanSoFar == 1,
%             fprintf('\t\tAnalyzing Channel 1 of %d',chanThisSess);
            numBack = length(sprintf('1 of %d',chanThisSess));
        else
            for ib = 1:numBack
%                 fprintf('\b');
            end
%             fprintf('%d of %d',chanSoFar,chanThisSess);
%             numBack = length(sprintf('%d of %d',chanSoFar,chanThisSess));
%             if chanSoFar == chanThisSess, fprintf('\n'); end
        end
        
        f = dir(sprintf('%s/%s/%s/DSP/%s/%s_%s_%s*',dataFold,monk{im},thisSess,thisChan,thisSess,thisChan,task));
        if length(f) == 0,
            %fprintf('\n\t\tNo "%s" file found for unit %s... Skipping to next unit\n',task,thisChan);
            outMat = cat(1,outMat,outRow);
            continue;
        elseif length(f) > 1,
            %fprintf('\n\t\tMultiple "%s" files found for unit %s... Please check directory structure\n',task,thisChan);
            %keyboard
        end
        fName = f(1).name;

        % Load unit data
        loadTic = tic;

        load(sprintf('%s/%s/%s/%s/%s/%s',dataFold,monk{im},thisSess,type,thisChan,fName));
        load(sprintf('%s/%s/%s/%s/%s/waves/%s_waves.mat',dataFold,monk{im},thisSess,type,thisChan,fName(1:(end-4))));
        % Next line can be used to time this load process, but clutters
        % the command line
%         fprintf('\t\tUnit data loaded at %s\n',printTiming(startTic));

        %% End learning section
        chanLoaded = chanLoaded + 1;
        
        % Check if the area string needs changed
        if sum(isletter(area)) == 0
            area = [' ',area];
        end
        area(strfind(area,'?')) = [];

        chanNum     = str2double(DSPname(4:end-1));
        
        % Get response property classification
        rf = getRF(spiketimes,Task);
        
        trialType = nan(size(spiketimes,1),1);
        trialType(Task.TargetLoc == rf) = 1;
        if isfield(Task,'DistLoc') && isfield(Task,'Singleton'), 
            trialType(Task.DistLoc == rf & Task.Singleton == 1) = 2;
            trialType(Task.DistLoc == rf & Task.Singleton == 0) = 3;
        end
        
%         fprintf('Getting stats at %s\n',printTiming(startTic));
        [pVals, pNames] = klGetTarget(spiketimes,Task);
        [sdf,sdfTimes]  = klSpkRatev2(spiketimes,'-q',1);
        
        targetLat       = klGetLatencyv4(spiketimes,Task,Task.TargetLoc,'-rasv',rf,{'StimOnset','SRT'},statSampRate,visualize);
        typeLat         = klGetLatencyv4(spiketimes,Task,trialType,'-rasv',rf,{'StimOnset','SRT'},statSampRate,visualize);
        
        if targetLat(1) > tstMax,
%             keyboard
            targetLat(1) = nan; 
        end
        if targetLat(2) < -tstMax || targetLat(2) > 30, 
%             keyboard
            targetLat(2) = nan; 
        end
        if typeLat(1) > tstMax, 
%             keyboard
            typeLat(1) = nan; 
        end
        if typeLat(2) < -tstMax || typeLat(2) > 30, 
%             keyboard
            typeLat(2) = nan; 
        end
        
%         if ir == 19,
%             keyboard
%         end
        
%         fprintf('Stats retrieved at %s\n',printTiming(startTic));
        
        outRow = {targetLat(1),targetLat(2),typeLat(1),typeLat(2),''};
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
%         fprintf('Row exported at %s\n',printTiming(startTic));
%         keyboard

    end
    xlswrite(outFile,outMat,monk{im},sprintf('%s5',num2abc(writeCol)));
    %keyboard
end


fprintf('***  Script completed in %s   ***\n',printTiming(startTic));

