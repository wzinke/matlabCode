function tdtPullSorts(file)

if ~exist('myDir','var'),
    myDir = 'Y:/Users/Kaleb/dataProcessed';
end
compType = 'euc';
minSpksTr = 5;
fresh = 0;

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

load(sprintf('%s/%s/sessSorts1.mat',myDir,file));

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
    
    sigSorts = find(sessSorts(ic).isSig);
    subAligned = chanSorts.alWaves(chanSorts.subInds,:);
    
    fprintf('\tExtrapolating Sort...');
    % Compare groups from subset to the whole set of waveforms
    switch compType,
        % Get means of the groups
        case 'corr', % Assign by maximum correlation
            mnSorts = nan(sessSorts(ic).k,size(chanSorts.alWaves,2));
            for i = 1:sessSorts(ic).k,
                mnSorts(i,:) = nanmean(subAligned(sessSorts(ic).idx==i,:),1);
            end
            allCorrs = corr(chanSorts.alWaves',mnSorts');
            [~,allIDX] = max(allCorrs,[],2);
        case 'euc', % Get minimum euclidean distance from centroid
            % Get transformation
            if ismember(sessSorts(ic).set,[1,2]),

                % If original data = score*coeff', then score =
                % data*inv(coeff')?
                allScores = chanSorts.alWaves*inv(chanSorts.pcaCoeffs');
                subScores = chanSorts.pca(:,1:2);
            else
                allScores = chanSorts.alWaves*chanSorts.lppEig;
                subScores = chanSorts.lpp(:,1:2);
            end

            % Get centroids
            mnSorts = nan(sessSorts(ic).k,size(subScores,2));
            for i = 1:sessSorts(ic).k,
                mnSorts(i,:) = nanmean(subScores(sessSorts(ic).idx==i,:),1);
            end
            sortDist = EuDist2(allScores(:,1:2),mnSorts);

            % Get minimum centroid distance
            [~,allIDX] = min(sortDist,[],2);
    end      
    
    for ib = 1:length(['Extrapolating Sort...']),
        fprintf('\b');
    end
    
    % Start sort loop
    if ~exist('trStarts','var'),
        trStarts = Task.trStarts;
        trEnds = Task.trEnds;
    end
    
    nTrs = length(trStarts);
    minSpks = minSpksTr*nTrs;
    
    goodUnits = 0;
    for is = 1:length(sigSorts),
        if ~fresh && exist(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',myDir,file,ic,is),'file'),
            fprintf('Unit already saved\n');
            continue;
        end
        clear sortTimes spikes
        fprintf('Counting unit %d spikes...',is);
        sortTimes = chanSorts.spkTimesAll(allIDX==sigSorts(is));
        spikes.allTimes = sortTimes;
        
        if length(sortTimes) < minSpks,
            continue;
        end
        goodUnits = goodUnits+1;
        
        % Figure out how big the spike matrix will be
        nSpks = nan(nTrs,1);
        for it = 1:nTrs,
            nSpks(it) = sum(sortTimes >= trStarts(it) & sortTimes <= trEnds(it));
        end
        maxSpk = max(nSpks(:));
        spikes.spiketimes = nan(nTrs,maxSpk);

        for ib = 1:length(sprintf('Counting unit %d spikes...',goodUnits)), fprintf('\b'); end
        fprintf('Placing unit %d spikes...',goodUnits);
        % Place spikes in the matrix
        for it = 1:nTrs,
            spikes.spiketimes(it,1:nSpks(it)) =  sortTimes(sortTimes >= trStarts(it) & sortTimes <= trEnds(it));
        end
        spikes.spiketimes = spikes.spiketimes-repmat(Task.AlignTimes,1,size(spikes.spiketimes,2));
        spikes.waves = chanSorts.alWaves(allIDX==sigSorts(is),:);
        for ib = 1:length(sprintf('Placing unit %d spikes...',goodUnits)), fprintf('\b'); end
        fprintf('Getting Spiking Statistics...');
        spikes.rateStruct = klRateStats(spikes.spiketimes);
        
        if ~exist(sprintf('%s/%s/Channel%d/Unit%d/',myDir,file,ic,goodUnits),'file'),
            mkdir(sprintf('%s/%s/Channel%d/Unit%d/',myDir,file,ic,goodUnits));
        end
        save(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',myDir,file,ic,goodUnits),'spikes');
        for ib = 1:length(sprintf('Getting Spiking Statistics...')), fprintf('\b'); end
        
    end
    fprintf('Channel %d completed\n',ic);
end
