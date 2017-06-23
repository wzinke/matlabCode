fresh = 0;

% Set up directories and certain constant filenames
xlFile = 'klTDTBookKeeping.xlsx';
procDir = 'Y:/Users/Kaleb/dataProcessed';
sessFile = 'tdtRecordingDoc.xlsx';

wvTimes = ((1:32)-9).*1000/24414;

% Get list of processed files
procFiles = dir(procDir);
procNames = {procFiles.name};
cutFolds = zeros(1,length(procNames));
for is = 1:length(procNames),
    cutFolds(is) = strcmp(procNames{is}(1),'.');
end
procNames(logical(cutFolds)) = [];

% Load in file
[excelNum,~,excelAll] = xlsread(xlFile,'All');
headerMat = excelAll(1:4,:);
if ~fresh,
    oldSess = unique(excelAll(5:end,1));
    procNames(ismember(procNames,oldSess)) = [];
    outMat = excelAll(5:end,:);
else
    % Initialize output matrix
    outMat = cell(0,size(excelAll,2));
end

% Start file loop
for is = 1:length(procNames),
    fprintf('Summarizing Session %s: ',procNames{is});
    
    % Load the excel infosheet for this session
    sessDate = procNames{is}((0:5)+find(ismember(procNames{is},'1'),1));
    [sessNum,~,sessAll] = xlsread(sessFile,sessDate);
    
    % Get some session info
    thisMonk = sessAll{strcmp(sessAll(:,1),'Monkey'),2};
    
    % Load "Task" for functional properties later
    load(sprintf('%s/%s/Behav.mat',procDir,procNames{is}));
    if ~isfield(Task,'trStarts'),
        Task.trStarts = trStarts;
        Task.trEnds = trEnds;
    end
    
    % Get probe info
    probeRow = find(strncmp(sessAll(:,1),'Probe',5),1);
    probeStarts = find(strncmp(sessAll(probeRow,:),'Probe',5));
    nProbes = length(probeStarts);
    areas = cell(1,nProbes); aps = cell(1,nProbes); mls = cell(1,nProbes); depths = nan(1,nProbes);
    for ip = 1:nProbes,
        areas{ip} = sessAll{strcmp(sessAll(:,probeStarts(ip)),'Area'),probeStarts(ip)+1};
        aps{ip} = sessAll{strcmp(sessAll(:,probeStarts(ip)),'AP/ML'),probeStarts(ip)+1};
        mls{ip} = sessAll{strcmp(sessAll(:,probeStarts(ip)),'AP/ML'),probeStarts(ip)+1};
        depths(ip) = sessAll{strcmp(sessAll(:,probeStarts(ip)),'Depth'),probeStarts(ip)+1};
    end
    
    % Get number of channels for loop
    chanDir = dir(sprintf('%s/%s/Channel*',procDir,procNames{is}));
    chanNames = {chanDir.name};
    nChans = length(chanNames);
    chans = nan(1,nChans);
    for i = 1:nChans,
        chans(i) = str2double(chanNames{i}(8:end));
    end
    chans = sort(chans);
    
    % Start channel loop
    for ic = 1:length(chans),
        fprintf('Channel %d of %d...',ic,length(chans));
        thisChan = chans(ic);
        thisProbe = floor((thisChan-1)/32)+1;
        
        % Get number of units here
        unitDir = dir(sprintf('%s/%s/Channel%d/Unit*',procDir,procNames{is},chans(ic)));
        unitNames = {unitDir.name};
        nUnits = length(unitNames);
        units = nan(1,nUnits);
        for i = 1:nUnits,
            units(i) = str2double(unitNames{i}(5:end));
        end
        units = sort(units);

        % Start unit loop
        for iu = 1:length(units),
            % Load in this unit's file
            load(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',procDir,procNames{is},chans(ic),units(iu)));
            outRow = cell(1,size(excelAll,2));
            
            %% Put things in an output row
            outRow{strcmp(excelAll(4,:),'Session')} = procNames{is};
            outRow{strcmp(excelAll(4,:),'Subject')} = thisMonk;
            outRow{strcmp(excelAll(4,:),'Channel')} = ic;
            outRow{strcmp(excelAll(4,:),'Unit')} = iu;
            
            outRow{strcmp(excelAll(4,:),'Area')} = areas{thisProbe};
            outRow{strcmp(excelAll(4,:),'AP')} = aps{thisProbe};
            outRow{strcmp(excelAll(4,:),'ML')} = mls{thisProbe};
            outRow{strcmp(excelAll(4,:),'Depth')} = depths(thisProbe);
            
            if ~isfield(spikes,'qualVect'),
                outRow{strcmp(excelAll(4,:),'SNR')} = klGetSNRv1(spikes.waves);
            else
                outRow{strcmp(excelAll(4,:),'SNR')} = spikes.qualVect(1);
                if length(spikes.qualVect) > 1,
                    outRow{strcmp(excelAll(4,:),'isoScore')} = spikes.qualVect(2);
                    outRow{strcmp(excelAll(4,:),'fnScore')} = spikes.qualVect(3);
                    outRow{strcmp(excelAll(4,:),'fpScore')} = spikes.qualVect(4);
                end
            end
            
            % Do waveform statistics
            [alWv,alTm] = klTroughAlignv5(spline(1:32,mean(spikes.waves,1),1:.1:32),spline(1:32,wvTimes,1:.1:32),0);
            [wvWidth, tfr] = klWvWidthv1(alWv,alTm);
            
            outRow{strcmp(excelAll(4,:),'wvWidth')} = wvWidth;
            outRow{strcmpi(excelAll(4,:),'TfR')} = tfr;
            
            % Get spiking statistics
            spkStarts = find(strcmp(excelAll(4,:),'mnRate'));
            spkTypes = fieldnames(spikes.rateStruct);
            for i = 1:length(spkTypes),
                spkVect = cell(1,10);
                spkVect{1} = spikes.rateStruct.(spkTypes{i}).mnRate;
                spkVect{2} = spikes.rateStruct.(spkTypes{i}).stdRate;
                spkVect{3} = spikes.rateStruct.(spkTypes{i}).fano;
                spkVect{4} = spikes.rateStruct.(spkTypes{i}).cv;
                spkVect{5} = spikes.rateStruct.(spkTypes{i}).cv2;
                spkVect{6} = spikes.rateStruct.(spkTypes{i}).lv;
                spkVect{7} = spikes.rateStruct.(spkTypes{i}).lvr;
                spkVect{8} = spikes.rateStruct.(spkTypes{i}).mnISI;
                spkVect{9} = spikes.rateStruct.(spkTypes{i}).stdISI;
                spkVect{10} = spikes.rateStruct.(spkTypes{i}).isiShort;
                
                outRow(spkStarts(i):(spkStarts(i)+length(spkVect)-1)) = spkVect;
            end
            
            % Now get functional properties
            if any(strcmpi(Task.TaskType,'MG')),
                [mgType, mgTune, mgLats] = klGetMG(spikes.spiketimes,Task);
                % Neuron type (vis, mov, vismov, none)
                outRow{strcmp(excelAll(4,:),'mgType')} = mgType;
                % Tuning properties
                outRow{strcmp(excelAll(4,:),'visDir')} = mgTune.vis.mu;
                outRow{find(strcmp(excelAll(4,:),'visDir'))+1} = [mgTune.vis.sig];
                outRow{find(strcmp(excelAll(4,:),'visDir'))+2} = mgTune.vis.amp;
                outRow{find(strcmp(excelAll(4,:),'visDir'))+3} = mgTune.vis.bl;
                outRow{strcmp(excelAll(4,:),'movDir')} = mgTune.mov.mu;
                outRow{find(strcmp(excelAll(4,:),'movDir'))+1} = mgTune.mov.sig;
                outRow{find(strcmp(excelAll(4,:),'movDir'))+2} = mgTune.mov.amp;
                outRow{find(strcmp(excelAll(4,:),'movDir'))+3} = mgTune.mov.bl;
                % Latencies
                outRow{find(strcmp(excelAll(4,:),'vRise'),1)} = mgLats.vis.rise;
                outRow{find(strcmp(excelAll(4,:),'vRise'),1)+1} = mgLats.vis.sig;
                outRow{find(strcmp(excelAll(4,:),'vRise'),1)+2} = mgLats.vis.poissLat;
                outRow{find(strcmp(excelAll(4,:),'vRise'),1)+3} = mgLats.vis.percPoiss;
                outRow{find(strcmp(excelAll(4,:),'mRise'),1)} = mgLats.mov.rise;
                outRow{find(strcmp(excelAll(4,:),'mRise'),1)+1} = mgLats.mov.sig;
                outRow{find(strcmp(excelAll(4,:),'mRise'),1)+2} = mgLats.mov.poissLat;
                outRow{find(strcmp(excelAll(4,:),'mRise'),1)+3} = mgLats.mov.percPoiss;
            end
            if any(ismember(Task.TaskType,{'Search','Cap','Capture'})),
                [capType,capTune,capLats] = klGetSearch(spikes.spiketimes,Task);
                % Neuron type (vis, mov, vismov, none) - Not the most
                % reliable but why not get it?
                outRow{strcmp(excelAll(4,:),'capType')} = capType;
                % Tuning properties
                outRow{strcmp(excelAll(4,:),'vRF')} = capTune.vRF;
                outRow{strcmp(excelAll(4,:),'mRF')} = capTune.mRF;
                % Latencies
                outRow{find(strcmp(excelAll(4,:),'vTargLoc'),1)}    = capLats.vis.targLocLat;
                outRow{find(strcmp(excelAll(4,:),'vTargLoc'),1)+1}  = capLats.vis.targTypeLat(1);
                outRow{find(strcmp(excelAll(4,:),'vTargLoc'),1)+2}  = capLats.mov.targTypeLat(1);
                outRow{find(strcmp(excelAll(4,:),'vTargLoc'),1)+3}  = capLats.vis.rise;
                outRow{find(strcmp(excelAll(4,:),'vTargLoc'),1)+4}  = capLats.vis.sig;
                outRow{find(strcmp(excelAll(4,:),'vTargLoc'),1)+5}  = capLats.vis.poissLat;
                outRow{find(strcmp(excelAll(4,:),'vTargLoc'),1)+6}  = capLats.vis.percPoiss;
                outRow{find(strcmp(excelAll(4,:),'mTargLoc'),1)}    = capLats.mov.targLocLat;
                outRow{find(strcmp(excelAll(4,:),'mTargLoc'),1)+1}  = capLats.vis.targTypeLat(2);
                outRow{find(strcmp(excelAll(4,:),'mTargLoc'),1)+2}  = capLats.mov.targTypeLat(2);
                outRow{find(strcmp(excelAll(4,:),'mTargLoc'),1)+3}  = capLats.mov.rise;
                outRow{find(strcmp(excelAll(4,:),'mTargLoc'),1)+4}  = capLats.mov.sig;
                outRow{find(strcmp(excelAll(4,:),'mTargLoc'),1)+5}  = capLats.mov.poissLat;
                outRow{find(strcmp(excelAll(4,:),'mTargLoc'),1)+6}  = capLats.mov.percPoiss;
            end
            outMat = cat(1,outMat,outRow);
            
        end
        for ib = 1:length(sprintf('Channel %d of %d...',ic,length(chans))), fprintf('\b'); end
        
    end
    fprintf('Done!\n');
end
xlswrite(xlFile,headerMat,'All','A1');
xlswrite(xlFile,outMat,'All','A5');

%% Now save individual monkeys
uMonks = unique(outMat(:,2));
for im = 1:length(uMonks),
    xlswrite(xlFile,headerMat,uMonks{im},'A1');
    xlswrite(xlFile,outMat(strcmp(outMat(:,2),uMonks{im}),:),uMonks{im},'A5');
end

            