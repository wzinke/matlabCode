function tdtSortSessionv2(inFile,varargin)

fresh = 0;
manThresh = 0;
nDims = 2;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-f','fresh'},
            fresh = varargin{varStrInd(iv)+1};
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

% Start channel loop
badSorts = [];
fprintf('\t\tFound %d channels\n',numChans);
for ic = 1:numChans,
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
%     try
        if manThresh,
            chanSorts = tdtSortChanv6a(chanSEVs.Wav1.data,'-t',chanThresh(ic),'-d',nDims);
        else
            chanSorts = tdtSortChanv6a(chanSEVs.Wav1.data,'-d',nDims);
        end

        % Save output in folder
        if ~exist(sprintf('%s/%s/Channel%d',procDir,inFile,ic),'file'),
            mkdir(sprintf('%s/%s/Channel%d',procDir,inFile,ic));
        end
        save(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',procDir,inFile,ic),'chanSorts','-v7.3');
%     catch
%         fprintf('Unable to sort channel %d\n',ic);
%         badSorts = cat(2,badSorts,ic);
%         
%     end
        
end
if ~isempty(badSorts),
    save(sprintf('%s/%s/Channel%d/badSorts.mat',procDir,inFile,ic),'badSorts');
end
        
        
        
        

    
 