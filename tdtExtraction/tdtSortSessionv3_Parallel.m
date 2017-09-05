function tdtSortSessionv3_Parallel(inFile,varargin)

fresh = 0;
manThresh = 0;
chanThresh = [];
nPar = 4;
nDims = 2;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-f','fresh'},
            fresh = varargin{varStrInd(iv)+1};
        case {'-t','thresh'},
            chanThresh = varargin{varStrInd(iv)+1};
        case {'-p','nPar'},
            nPar = varargin{varStrInd(iv)+1};
        case {'-d','nDims'},
            nDims = varargin{varStrInd(iv)+1}; 
    end
end

% Set directories
rawDir = 'Y:/Users/Kaleb/dataRaw';
procDir = 'Y:/Users/Kaleb/dataProcessed';

% Get processed session names
fprintf('\tSorting file %s:\n',inFile);

% Set up processed folder
if ~exist(sprintf('%s/%s',procDir,inFile),'file'),    
    mkdir(sprintf('%s/%s',procDir,inFile));
end

% Get number of channels to deal with
numChans = length(dir(sprintf('%s/%s/*.sev',rawDir,inFile)))/2;

if exist(sprintf('%s/%s/manChanThreshes.mat',procDir,inFile),'file'),
    load(sprintf('%s/%s/manChanThreshes.mat',procDir,inFile));
    manThresh = 1;
end

% Open up a parallel pool, if necessary
x=tic;
if isempty(gcp('nocreate')),
    parPool = parpool(nPar);
end
fprintf('\n\nPool opened in %s\n\n',printTiming(x));


% Start channel loop
fprintf('\t\tFound %d channels\n',numChans);
badSorts = nan(1,numChans);
% [allSorts(1:numChans).neg] = deal(struct('reConstruct',[],'subInds',[],'idx',[],'aggMap',[],'allWaves',[],'k',[],'isSig',[],'subScores',[],'dimRed',[],'allTimes',[],'randStruct',[]));
% [allSorts(1:numChans).pos] = deal(struct('reConstruct',[],'subInds',[],'idx',[],'aggMap',[],'allWaves',[],'k',[],'isSig',[],'subScores',[],'dimRed',[],'allTimes',[],'randStruct',[]));

parfor ic = 1:numChans,
%     clear chanSEVs chanSorts
    if ~fresh && exist(sprintf('%s/%s/Channel%d/autoSort_noAudit.mat',procDir,inFile,ic),'file'),
        fprintf('\t\t\tChannel %d is already sorted... Moving on\n',ic);
        continue;
    else
        fprintf('\t\t\tReading channel %d: ',ic);
    end
    
    % Read in data
    chanSEVs = SEV2mat_kl(sprintf('%s/%s',rawDir,inFile),'CHANNEL',ic,'VERBOSE',0,'EVENTNAME','Wav1');
    
    while ~any(chanSEVs.Wav1.data > 1),
        chanSEVs.Wav1.data = chanSEVs.Wav1.data.*1000;
        chanThresh(ic) = chanThresh(ic).*1000;
    end
    
    % Do the sorting
    try
        if manThresh,
            chanSorts = tdtSortChanv6a(chanSEVs.Wav1.data,'-t',chanThresh(ic),'-d',nDims);
        else
            chanSorts = tdtSortChanv6a(chanSEVs.Wav1.data,'-d',nDims);
        end

        % Save output in folder
        if ~exist(sprintf('%s/%s/Channel%d',procDir,inFile,ic),'file'),
            mkdir(sprintf('%s/%s/Channel%d',procDir,inFile,ic));
        end
        saveStr = sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',procDir,inFile,ic);
        while exist(saveStr,'file'),
            saveStr = [saveStr,'_'];
        end
        klParSavev1(saveStr,chanSorts,'chanSorts');
    catch
        fprintf('Unable to sort channel %d\n',ic);
        badSorts(ic) = 1;
    end
        
end
if ~isempty(badSorts),
    save(sprintf('%s/%s/Channel%d/badSorts.mat',procDir,inFile),'badSorts');
end
        
delete(parPool);        


 