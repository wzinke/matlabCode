% % tdtExtractMaster

fprintf('*********************************************\n');
fprintf('****** BEGINNING MASTER TDT EXTRACTION ******\n');
fprintf('*********************************************\n\n\n');

% Set constants
lfpWind = [-500, 6500];
lfpFreq = 1000/1017; % 1000 converts to ms
spkFreq = 1000/24414;
fresh = 0;

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
procNames(ismember(procNames,{'.','..'})) = [];

% Check for unprocessed files
if ~fresh,
    for ir = 1:length(rawNames),
        isExtracted(ir) = tdtIsExtracted(rawNames{ir});
    end
    toDo = rawNames(~isExtracted);
    fprintf('**** FOUND %d UNPROCESSED FILES ****\n',length(toDo));
else
    toDo = rawNames;
    fprintf('*********  FOUND %d  FILES *********\n',length(toDo));
end

% Load up events
events=TEMPO_EV_cosman_rig028;

% Start session loop
for is = 1:length(toDo),
    thisFile = toDo{is};
    fprintf('\tProcessing file %s:\n',thisFile);
    success = 0;
    
    % Set up processed folder
    mkdir(sprintf('%s/%s',procDir,thisFile));
    
    % First, get events and task stuff
    fprintf('\t\tBehavior extraction: ',thisFile);
    tdtExtractSessv4(thisFile,events);
    
    success = 1;
end
            
        