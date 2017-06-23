function tdtAskThresh(inFile)

if ~exist('rawDir','var') || isempty(rawDir),
    rawDir = 'Y:/Users/Kaleb/dataRaw';
end
if ~exist('procDir','var') || isempty(procDir),
    procDir = 'Y:/Users/Kaleb/dataProcessed';
end

spkFreq = 1000/24414;


% Now, channel by channel, get spikes and LFPs
numChans = length(dir(sprintf('%s/%s/*.sev',rawDir,inFile)))/2;
fprintf('\t\tNeural Data: Found %d channels\n',numChans);

mySpks = nan(numChans,ceil(5000/spkFreq));
for ic = 1:numChans,
    printStr = sprintf('\t\t\tFetching cell %d of %d...',ic,numChans);
    fprintf(printStr);
    clear chanSEVs
    
    chanSEVs = SEV2mat_kl(sprintf('%s/%s/',rawDir,inFile),'CHANNEL',ic,'VERBOSE',0);
    spkTimesRaw = 0:spkFreq:(spkFreq*(size(chanSEVs.Wav1.data,2)-1)); % 1000x multiplier converts to ms

    randStart = randi((length(spkTimesRaw)-ceil((5000/spkFreq))));
    randEnd = randStart + ceil(5000/spkFreq) - 1;
    
    mySpks(ic,:) = chanSEVs.Wav1.data(randStart:randEnd);
    for ib = 1:length(printStr),
        fprintf('\b');
    end
end
fprintf('\n');
for ic = 1:numChans,
    printStr = sprintf('\t\t\tShowing cell %d of %d...',ic,numChans);
    fprintf(printStr);

    figure(ic);
    plot(mySpks(ic,:)');
    
    picked = 0;
    while ~picked,
        try
            chanThresh(ic) = input('Please enter the threshold to use for this channel: ');
            picked = 1;
        end
    end
    close(ic);
    
    for ib = 1:(length(printStr)+length('Please enter the threshold to use for htis channel: ')),
        fprintf('\b');
    end
end
if ~exist(sprintf('%s/%s',procDir,inFile),'file'),
    mkdir(sprintf('%s/%s',procDir,inFile));
end
save(sprintf('%s/%s/manChanThreshes.mat',procDir,inFile),'chanThresh');
    