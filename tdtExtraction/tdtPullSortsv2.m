function tdtPullSortsv2(file)

if ~exist('myDir','var'),
    myDir = 'Y:/Users/Kaleb/dataProcessed';
end
compType = 'euc';
minSpksTr = 5;
fresh = 0;
myPols = {'neg','pos'};
maxSpkTm = 5000;

% Get number of channels for loop
chanDir = dir(sprintf('%s/%s/Channel*',myDir,file));
chanNames = {chanDir.name};
nChans = length(chanNames);
chans = nan(1,nChans);
for i = 1:nChans,
    chans(i) = str2num(chanNames{i}(8:end));
end
chans = sort(chans);

% Load task if necessary
if nChans > 0,
    load(sprintf('%s/%s/Behav.mat',myDir,file));
end

% load(sprintf('%s/%s/sessSorts1.mat',myDir,file));

% Start channel loop
for ic = 1:nChans,
    chanFold = sprintf('%s/%s/Channel%d',myDir,file,chans(ic));
    
    if exist(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',myDir,file,ic),'file'),
        fprintf('Loading Channel %d...',ic);
        load(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',myDir,file,chans(ic)));
        fprintf('Done!\n');
    else
        fprintf('Channel %d not yet sorted\n',ic);
        continue
    end
    
    if isempty(chanSorts.neg.allWaves),
        fprintf('No threshold crossings on channel %d... Moving on\n',ic);
        continue
    end
    
    % Start sort loop
    if ~exist('trStarts','var'),
        trStarts = Task.trStarts;
        trEnds = Task.trEnds;
    end
    
    nTrs = length(trStarts);
    minSpks = minSpksTr*nTrs;
    
    if ~isfield(chanSorts.neg,'idx') || ~isfield(chanSorts.pos,'idx'),
        fprintf('*** Channel %d not sorted... Moving on ***');
        continue
    end
    
    nUnits = chanSorts.neg.k + chanSorts.pos.k;
    unitPol = [ones(1,chanSorts.neg.k), ones(1,chanSorts.pos.k).*2];
    unitID = [1:chanSorts.neg.k,1:chanSorts.pos.k];
    for is = 1:nUnits,
        if ~fresh && exist(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',myDir,file,ic,is),'file'),
            fprintf('Unit already saved\n');
            continue;
        end
        clear sortTimes spikes
        fprintf('Counting unit %d spikes...',is);
        sortTimes = chanSorts.(myPols{unitPol(is)}).allTimes(chanSorts.(myPols{unitPol(is)}).allIDX==unitID(is));
        spikes.allTimes = sortTimes;
        
        if length(sortTimes) < minSpks,
            continue;
        end
%         goodUnits = goodUnits+1;
        
        % Figure out how big the spike matrix will be
        nSpks = nan(nTrs,1);
        for it = 1:nTrs,
            nSpks(it) = sum(sortTimes >= trStarts(it) & sortTimes <= trEnds(it));
        end
        maxSpk = max(nSpks(:));
        spikes.spiketimes = nan(nTrs,maxSpk);

        for ib = 1:length(sprintf('Counting unit %d spikes...',is)), fprintf('\b'); end
        fprintf('Placing unit %d spikes...',is);
        % Place spikes in the matrix
        for it = 1:nTrs,
            spikes.spiketimes(it,1:nSpks(it)) =  sortTimes(sortTimes >= trStarts(it) & sortTimes <= trEnds(it));
        end
        spikes.spiketimes = spikes.spiketimes-repmat(Task.AlignTimes,1,size(spikes.spiketimes,2));
        spikes.spiketimes(spikes.spiketimes >= maxSpkTm) = nan;
        cutCols = zeros(1,size(spikes.spiketimes,2));
        for iCol = 1:size(spikes.spiketimes,2)
            cutCols(iCol) = sum(isnan(spikes.spiketimes(:,iCol)))==size(spikes.spiketimes,1);
        end
        spikes.spiketimes(:,logical(cutCols)) = [];
        
        spikes.waves = chanSorts.(myPols{unitPol(is)}).allWaves(chanSorts.(myPols{unitPol(is)}).allIDX==unitID(is),:);
        for ib = 1:length(sprintf('Placing unit %d spikes...',is)), fprintf('\b'); end
        fprintf('Getting Spiking Statistics...');
        spikes.rateStruct = klRateStats(spikes.spiketimes);
        
        qualVect(1) = klGetSNRv1(chanSorts.(myPols{unitPol(is)}).allWaves(chanSorts.(myPols{unitPol(is)}).allIDX==unitID(is),:));
%         [qualVect(2),simMat] = klGetISv2(chanSorts.(myPols{unitPol(is)}).normObs(chanSorts.(myPols{unitPol(is)}).subInds(chanSorts.(myPols{unitPol(is)}).idx(:,chanSorts.(myPols{unitPol(is)}).k)==unitID(is)),:),chanSorts.(myPols{unitPol(is)}).normObs(chanSorts.(myPols{unitPol(is)}).subInds(chanSorts.(myPols{unitPol(is)}).idx(:,chanSorts.(myPols{unitPol(is)}).k)~=unitID(is)),:));
%         qualVect(3) = klGetFNv2(chanSorts.(myPols{unitPol(is)}).normObs(chanSorts.(myPols{unitPol(is)}).subInds(chanSorts.(myPols{unitPol(is)}).idx(:,chanSorts.(myPols{unitPol(is)}).k)==unitID(is)),:),chanSorts.(myPols{unitPol(is)}).normObs(chanSorts.(myPols{unitPol(is)}).subInds(chanSorts.(myPols{unitPol(is)}).idx(:,chanSorts.(myPols{unitPol(is)}).k)~=unitID(is)),:),'sim',simMat);
%         qualVect(4) = klGetFPv2(chanSorts.(myPols{unitPol(is)}).normObs(chanSorts.(myPols{unitPol(is)}).subInds(chanSorts.(myPols{unitPol(is)}).idx(:,chanSorts.(myPols{unitPol(is)}).k)==unitID(is)),:),chanSorts.(myPols{unitPol(is)}).normObs(chanSorts.(myPols{unitPol(is)}).subInds(chanSorts.(myPols{unitPol(is)}).idx(:,chanSorts.(myPols{unitPol(is)}).k)~=unitID(is)),:),'sim',simMat);
        spikes.qualVect = qualVect;
        
        if ~exist(sprintf('%s/%s/Channel%d/Unit%d/',myDir,file,ic,is),'file'),
            mkdir(sprintf('%s/%s/Channel%d/Unit%d/',myDir,file,ic,is));
        end
        save(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',myDir,file,ic,is),'spikes');
        for ib = 1:length(sprintf('Getting Spiking Statistics...')), fprintf('\b'); end
        
    end
    fprintf('Channel %d completed\n',ic);
end
