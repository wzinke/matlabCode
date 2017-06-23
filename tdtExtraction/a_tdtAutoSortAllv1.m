% % % TDT Master Extract, pre-sort
globalTic = tic;
fprintf('*********************************************\n');
fprintf('****** BEGINNING MASTER TDT EXTRACTION ******\n');
fprintf('*********************************************\n\n\n');

% Set constants
lfpWind = [-500, 2500];
lfpFreq = 1000/1017; % 1000 converts to ms
spkFreq = 1000/24414;
fresh = 0;
doPar = 1;
manual = 0;
defThresh = 35;
nPar = 2;
nDims = 2;

% Set directories
rawDir = 'Y:/Users/Kaleb/dataRaw';
procDir = 'Y:/Users/Kaleb/dataProcessed';
sessFile = 'C:/Users/Kaleb/Documents/MATLAB/tdtRecordingDoc.xlsx';

% Get raw session names
rawSess = dir(rawDir);
rawNames = {rawSess.name};
rawNames(ismember(rawNames,{'.','..','.AppleDouble'})) = [];

% Get processed session names
procSess = dir(procDir);
procNames = {procSess.name};
procNames(ismember(procNames,{'.','..','.AppleDouble'})) = [];

% Check for unprocessed files
if ~fresh,
%     for ir = 1:length(rawNames),
%         isExtracted(ir) = tdtIsExtracted(rawNames{ir});
%     end
%     toDo = rawNames(~isExtracted);
    toDo = rawNames(~ismember(rawNames,procNames));
    fprintf('**** FOUND %d UNPROCESSED FILES ****\n',length(toDo));
else
    toDo = rawNames;
    fprintf('*********  FOUND %d  FILES *********\n',length(toDo));
end

% Load up events
events=TEMPO_EV_cosman_rig028;

% Start session loop

% Get thresholds for sorting
for is = 1:length(toDo),
    if ~exist(sprintf('%s/%s/manChanThreshes.mat',procDir,toDo{is}),'file'),
        if manual,
            tdtAskThresh(toDo{is});
        else
            sessDate = toDo{is}((0:5)+find(ismember(toDo{is},'1'),1));
            [sessNum,~,sessAll] = xlsread(sessFile,sessDate);
            probStart = find(strcmp(sessAll(9,:),'P2P'));
            chanThresh = [];
            for ip = 1:length(probStart),
                for ic = 1:32,
                    if ischar(sessAll{9+ic,probStart(ip)}),
                        slashInd = strfind(sessAll{9+ic,probStart(ip)},'/');
                        denom = str2double(sessAll{9+ic,probStart(ip)}((slashInd+1):end));
                        chanThresh = cat(2,chanThresh,ceil(denom/2));
                    else
                        chanThresh = cat(2,chanThresh,defThresh);
                    end
                end
            end
            if ~exist(sprintf('%s/%s',procDir,toDo{is}),'file'),
                mkdir(sprintf('%s/%s',procDir,toDo{is}));
            end
            save(sprintf('%s/%s/manChanThreshes.mat',procDir,toDo{is}),'chanThresh');
        end
    end
end

% Now extract channels
for is = 1:length(toDo),
    clearvars -except events toDo procDir rawDir is globalTic doPar nPar
    thisFile = toDo{is};
    fprintf('\tProcessing file %s:\n',thisFile);
    success = 0;
    
    % Set up processed folder
    mkdir(sprintf('%s/%s',procDir,thisFile));
    
    % First, get events and task stuff
    fprintf('\t\tBehavior extraction: ',thisFile);
    Task = tdtGetTaskv4(sprintf('%s/%s',rawDir,thisFile),events);
    
    % Save behavior only
    fprintf('Saving...')
    save(sprintf('%s/%s/Behav.mat',procDir,thisFile),'Task');
    for ib = 1:length(['Saving...']), fprintf('\b'); end
    fprintf('Done!\n');
    
    % Now, make a folder for each channel and get MUA/LFPs
    numChans = length(dir(sprintf('%s/%s/*.sev',rawDir,thisFile)))/2;
    fprintf('\t\tNeural Data: Found %d channels\n',numChans);
    for ic = 1:numChans,
%         clear spikes LFP
%         [spikes, LFP] = tdtGetChanv2(sprintf('%s/%s',rawDir,thisFile),ic,Task);
%         
        mkdir(sprintf('%s/%s/Channel%d',procDir,thisFile,ic));
%         save(sprintf('%s/%s/Channel%d/Spikes.mat',procDir,thisFile,ic),'Task','spikes');
%         save(sprintf('%s/%s/Channel%d/LFPs.mat',procDir,thisFile,ic),'Task','LFP');
    end
    % Now we'll sort this session...
    if exist(sprintf('%s/%s/manChanThreshes.mat',procDir,thisFile),'file'),
        load(sprintf('%s/%s/manChanThreshes.mat',procDir,thisFile));
    end
    if doPar,
        if exist('chanThresh','var'),
            tdtSortSessionv3_Parallel(thisFile,'-t',chanThresh,'-p',nPar);
        else
            tdtSortSessionv3_Parallel(thisFile,'-p',nPar);
        end
    else
        if exist('chanThresh','var'),
            tdtSortSessionv2(thisFile,'-t',chanThresh);
        else
            tdtSortSessionv2(thisFile);
        end
    end
    
    % Take results and save them in unit by unit folders
    tdtPullSorts_Agglomv2(thisFile);
end

% End: Now move on to actually sorting them and wrapping up in step 2
    
fprintf('******     STEP 1 COMPLETE IN %s      ********\n',printTiming(globalTic));
    
    