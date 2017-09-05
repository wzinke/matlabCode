fresh = 0;

% Set directories
rawDir = 'Y:/Users/Kaleb/dataRaw';
procDir = 'Y:/Users/Kaleb/dataProcessed';

% Get step-1 processed session names
procSess = dir(procDir);
procNames = {procSess.name};
procNames(ismember(procNames,{'.','..','.AppleDouble'})) = [];

% Check for unprocessed files
if ~fresh,
    isSorted = zeros(1,length(procNames));
    for ip = 1:length(procNames),
        isSorted(ip) = any(exist(sprintf('%s/%s/sessSorts1.mat',procDir,procNames{ip}),'file'));
    end
    toDo = procNames(~isSorted);
    fprintf('**** FOUND %d UNSORTED FILES ****\n',length(toDo));
else
    toDo = procNames;
    fprintf('*********  FOUND %d  FILES *********\n',length(toDo));
end

% Now start the session by session sorting loop
for is = 1:length(toDo),
    clearvars -except is toDo procDir rawDir
    thisFile = toDo{is};
    
    % Open and view potential sorts, select appropriate ones
    tdtViewSorts(thisFile);
    
    % Take results and save them in unit by unit folders
    tdtPullSorts(thisFile);
    
    % Now get sort qualities
    tdtGetQuals(thisFile);
end