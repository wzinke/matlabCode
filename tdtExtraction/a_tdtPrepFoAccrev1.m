% % % TDT Master Extract, pre-sort
globalTic = tic;
fprintf('*********************************************\n');
fprintf('****** BEGINNING MASTER TDT EXTRACTION ******\n');
fprintf('*********************************************\n\n\n');

% Set SLURM params
nCores =  1;
nTasks = 1;
memPerCore = 64;
memUnits = 'G';
outFilePref = 'myJob_';
nNodes = 1;
timeStr = '04:00:00';
slurmShellName = 'submitAllSlurm.txt';
remDir = './data';

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
maxWvs = 30000;
extrapMax = 10000;
minRefract = .6;
doSimple = 0;
transferAccre = 1;
submitAccre = 1;
mailType = 'FAIL';

% Set directories
rawDir = 'Y:/Users/Kaleb/dataRaw';
procDir = 'Y:/Users/Kaleb/dataProcessed';
sessFile = 'C:/Users/Kaleb/Documents/MATLAB/tdtExtraction/tdtRecordingDoc.xlsx';

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

% Start session loop and hide username/password a little bit...
userOpts.name = 'loweka';
userOpts.pass = 'V@ndy20!3k';

% Get thresholds for sorting
for is = 1:length(toDo),
    if ~exist(sprintf('%s/%s/manChanThreshes.mat',procDir,toDo{is}),'file'),
        if manual,
            tdtAskThresh(toDo{is});
        else
            sessDate = toDo{is}((0:5)+find(ismember(toDo{is},'1'),1));
            [sessNum,~,sessAll] = xlsread(sessFile,sessDate);
            probStart = find(strcmp(sessAll(11,:),'P2P'));
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
for is = 1%:length(toDo),
    clearvars -except mailType slurmShellName maxWvs nDims nCores memPerCore...
        memUnits outFilePref userOpts transferAccre submitAccre...
        timeStr minRefract spkFreq doSimple events toDo procDir rawDir...
        is globalTic doPar nPar remDir extrapMax
    thisFile = toDo{is};
    fprintf('\tProcessing file %s:\n',thisFile);
    success = 0;
    
    % Set up processed folder
    mkdir(sprintf('%s/%s',procDir,thisFile));
    
    % First, get events and task stuff
    fprintf('\t\tBehavior extraction: ',thisFile);
    try
        Task = tdtGetTaskv4(sprintf('%s/%s',rawDir,thisFile),events);
        
        % Get eye position data
        eyes = tdtGetEyes(sprintf('%s/%s',rawDir,thisFile),Task);
        
        % Save behavior only
        fprintf('Saving...')
        save(sprintf('%s/%s/Behav.mat',procDir,thisFile),'Task');
        save(sprintf('%s/%s/Eyes.mat',procDir,thisFile),'eyes');
        for ib = 1:length(['Saving...']), fprintf('\b'); end
        fprintf('Done!\n');
    catch
        fprintf('Unable to extract Task information... skipping this section\n');
    end
    
    % Get number of channels
    numChans = length(dir(sprintf('%s/%s/*.sev',rawDir,thisFile)))/2;
    fprintf('\t\tNeural Data: Found %d channels\n',numChans);
    
    if exist(sprintf('%s/%s/manChanThreshes.mat',procDir,thisFile),'file'),
        load(sprintf('%s/%s/manChanThreshes.mat',procDir,thisFile));
    end
    if ~exist('chanThresh','var'), chanThresh = ones(1,numChans).*defThresh; end
    
    % Now, make a folder for each channel and get MUA/LFPs
    notSaved = zeros(1,numChans);
	shellSubmit = fopen(sprintf('%s/%s/%s',procDir,thisFile,slurmShellName),'w+');
    for ic = 1:numChans,
        fprintf('\t\t\tWorking on channel %d: ',ic);
%         clear spikes LFP
        mkdir(sprintf('%s/%s/Channel%d',procDir,thisFile,ic));
        if doSimple,
            [spikes, LFP] = tdtGetChanv2(sprintf('%s/%s',rawDir,thisFile),ic,Task);
%         
            save(sprintf('%s/%s/Channel%d/Spikes.mat',procDir,thisFile,ic),'Task','spikes');
            save(sprintf('%s/%s/Channel%d/LFPs.mat',procDir,thisFile,ic),'Task','LFP');
            spikes = []; LFP = [];
        end
        printStr = 'Loading channel data...';
        fprintf('%s',printStr);
        chanData = SEV2mat_kl(sprintf('%s/%s',rawDir,thisFile),'VERBOSE',0,'CHANNEL',ic,'EVENTNAME','Wav1');
        while ~any(chanData.Wav1.data > 1),
            chanData.Wav1.data = chanData.Wav1.data.*1000;
        end
        chanTimes = single(0:spkFreq:(spkFreq*(size(chanData.Wav1.data,2)-1)));
        
        % Get threshold crossings
        for ib = 1:length(printStr), fprintf('\b'); end
        printStr = 'Getting threshold crossings...';
        fprintf('%s',printStr);
        [spkTimes,spkWaves,spkPols] = klThreshCrossv6(chanData.Wav1.data,'times',chanTimes,'-m',minRefract,'-t',chanThresh(ic));
        
        % clear chanData
        chanData = []; chanTimes = [];
        
        % Store waveforms and spike times
        for ib = 1:length(printStr), fprintf('\b'); end
        printStr = 'Storing crossings in structure...';
        fprintf('%s',printStr);
        chanSorts.neg.allWaves = single(spkWaves(spkPols==2,:));
        chanSorts.neg.allTimes = single(spkTimes(spkPols==2));
        chanSorts.pos.allWaves = single(spkWaves(spkPols==1,:));
        chanSorts.pos.allTimes = single(spkTimes(spkPols==1));
        spkTimes = []; spkWaves = []; spkPols = [];
        
        % Save it
        for ib = 1:length(printStr), fprintf('\b'); end
        printStr = 'Saving structure...';
        fprintf('%s',printStr);
        % This sometimes makes an error saving as v7.3, so let's try to
        % handle that
        try
            save(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',procDir,thisFile,ic),'chanSorts','-v7.3');
        catch saveErr
            if strcmpi(saveErr.identifier,'MATLAB:save:cantWriteFile'),
                save(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',procDir,thisFile,ic),'chanSorts');
            else
                notSaved(ic) = 1;
            end
        end
        chanSorts = [];
        
        % Make SLURM file
        [slurmPath,thisSlurm] = klMakeSlurmv3('klChanSortv1.m','-f',...
			sprintf('''%s'',%d,''-w'',%d,''-d'',%d,''-r'',''%s'',''-extrap'',%d',thisFile,ic,maxWvs,nDims,remDir,extrapMax),...
            '-t',timeStr,'-s',sprintf('%s_Ch%d.slurm',thisFile,ic),'sDir',sprintf('%s/%s',procDir,thisFile),...
            '-c',nCores,'-mpc',memPerCore,'-mu',memUnits,'-m',memPerCore,'-out',[outFilePref,thisFile,'_Ch',num2str(ic),'.out'],'-mt',mailType);
        
		fprintf(shellSubmit,sprintf('sbatch %s/%s/%s',remDir,thisFile,thisSlurm));
		fprintf(shellSubmit,'\n');
        
        % Confirm that it's finished
        for ib = 1:length(printStr), fprintf('\b'); end
        printStr = 'Done!';
        fprintf('%s',printStr); fprintf('\n');
    end
    fclose(shellSubmit);
	
    % Now, send this session's folder in dataProcessed to ACCRE
    if transferAccre,
        success = klSendFoldToAccre(sprintf('%s/%s',procDir,thisFile),thisFile,userOpts,'-rdir',remDir,'-ldir','./tempTars');
    end
    if submitAccre && success,
        sshStruct = ssh2_config('login.accre.vanderbilt.edu',userOpts.name,userOpts.pass);
        sshStruct = ssh2_command(sshStruct,sprintf('sh %s/%s/%s',remDir,thisFile,slurmShellName));
        sshStruct = ssh2_close(sshStruct);
    end
end
    
fprintf('******     STEP 1 COMPLETE IN %s      ********\n',printTiming(globalTic));
    
    