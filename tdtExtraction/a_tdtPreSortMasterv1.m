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
doPar = 0;

% Set directories
rawDir = 'Y:/Users/Kaleb/dataRaw';
procDir = 'Y:/Users/Kaleb/dataProcessed';

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
% for is = 1:length(toDo),
%     tdtAskThresh(toDo{is});
% end

% Now extract channels
for is = 1:length(toDo),
    clearvars -except events toDo procDir rawDir is globalTic doPar
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
            tdtSortSessionv2_Parallel(thisFile,'-t',chanThresh);
        else
            tdtSortSessionv2_Parallel(thisFile);
        end
    else
        if exist('chanThresh','var'),
            tdtSortSessionv2(thisFile,'-t',chanThresh);
        else
            tdtSortSessionv2(thisFile);
        end
    end
    
    % Take results and save them in unit by unit folders
    tdtPullSorts_Agglom(thisFile);
    
    % Now get sort qualities
    tdtGetQuals(thisFile);
    
end

% End: Now move on to actually sorting them and wrapping up in step 2
    
fprintf('******     STEP 1 COMPLETE IN %s      ********\n',printTiming(globalTic));
    
    