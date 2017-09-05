function tdtPullSorts_Agglom(file)

if ~exist('myDir','var'),
    myDir = 'Y:/Users/Kaleb/dataProcessed';
end
compType = 'euc';
minSpksTr = 5;
maxSpkTm = 10000;
fresh = 0;
minCuts = 100;
myPols = {'neg','pos'};
sdStart = 3;

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
%     fprintf('Showing channel %d...\n',chans(ic));
    chanFold = sprintf('%s/%s/Channel%d',myDir,file,chans(ic));
%     close all;
%     clear sortAudit chanSorts subWaves
    
    if exist(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',myDir,file,ic),'file'),
        fprintf('Loading Channel %d...',ic);
        load(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',myDir,file,chans(ic)));
        fprintf('Done!\n');
    else
        fprintf('Channel %d not yet sorted\n',ic);
        continue
    end
    
    allIDX = [];
    allTimes = [];
    allPols = [];
    allSigs = [];
    sigPols = [];
    for ip = 1:length(myPols),
        sigSorts = find(chanSorts.(myPols{ip}).isSig);
        subAligned = chanSorts.(myPols{ip}).allWaves(chanSorts.(myPols{ip}).subInds,:);

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
                if strcmpi(chanSorts.(myPols{ip}).dimRed,'pca'),

                    % If original data = score*coeff', then score =
                    % data*inv(coeff')?
                    allScores = chanSorts.(myPols{ip}).allWaves*inv(chanSorts.(myPols{ip}).reConstruct');
                    subScores = chanSorts.(myPols{ip}).subScores(:,1:2);
                else
                    allScores = chanSorts.(myPols{ip}).allWaves*chanSorts.(myPols{ip}).reConstruct;
                    subScores = chanSorts.(myPols{ip}).subScores(:,1:2);
                end
                
                if length(unique(chanSorts.(myPols{ip}).idx)) > 1,
                    keyboard
                end
                
                % Get centroids
                mnSorts = nan(chanSorts.(myPols{ip}).k,size(subScores,2));
                for i = 1:chanSorts.(myPols{ip}).k,
                    mnSorts(i,:) = nanmean(subScores(chanSorts.(myPols{ip}).idx==i,:),1);
                end
                if any(isnan(chanSorts.(myPols{ip}).idx)),
                    mnSorts = cat(1,mnSorts,nanmean(subScores(isnan(chanSorts.(myPols{ip}).idx),:),1));
                end
                sortDist = EuDist2(allScores(:,1:2),mnSorts);

                % Get minimum centroid distance
                [~,polIDX] = min(sortDist,[],2);
                if length(unique(polIDX)) > 1,
                    polIDX(polIDX==max(polIDX))= nan;
                    for is = 1:length(unique(polIDX(~isnan(polIDX)))),
                        nCuts = inf;
                        thisSD = sdStart;
                        while nCuts > minCuts,
                            thisZ = (allScores(:,1:2)-repmat(nanmean(allScores(polIDX==is,1:2),1),size(allScores,1),1))./repmat(nanstd(allScores(polIDX==is,1:2),[],1),size(allScores,1),1);
                            zRad = sqrt(sum(thisZ.^2,2));
                            nCuts = sum(polIDX==is & zRad > thisSD);
                            polIDX(polIDX==is & zRad > thisSD) = nan;
                            thisSD = thisSD+1;
                        end
                    end
%                     keyboard
                end
                
                polTimes = chanSorts.(myPols{ip}).allTimes';
                polID = ones(length(polTimes),1).*ip;
                
        end      
        allIDX{ip} = polIDX;
        allTimes = cat(1,allTimes,polTimes);
        allPols = cat(1,allPols,polID);
        allSigs = cat(2,allSigs,sigSorts);
        sigPols = cat(2,sigPols,ones(1,length(sigSorts)).*ip);
%         keyboard
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
    for is = 1:length(allSigs),
        if ~fresh && exist(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',myDir,file,ic,is),'file'),
            fprintf('Unit already saved\n');
            continue;
        end
        clear sortTimes spikes
        fprintf('Counting unit %d spikes...',is);
        sortTimes = chanSorts.(myPols{sigPols(is)}).allTimes(allIDX{sigPols(is)}==allSigs(is));
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
        spikes.spiketimes(spikes.spiketimes >= maxSpkTm) = nan;
        spikes.waves = chanSorts.(myPols{sigPols(is)}).allWaves(allIDX{sigPols(is)}==allSigs(is),:);
        for ib = 1:length(sprintf('Placing unit %d spikes...',goodUnits)), fprintf('\b'); end
        fprintf('Getting Spiking Statistics...');
        spikes.rateStruct = klRateStats(spikes.spiketimes);

        % Now get sort qualities
        qualVect(1) = klGetSNRv1(chanSorts.(myPols{sigPols(is)}).allWaves(allIDX{sigPols(is)}==allSigs(is),:));
        [qualVect(2),simMat] = klGetISv2(chanSorts.(myPols{sigPols(is)}).subScores(chanSorts.(myPols{sigPols(is)}).idx==allSigs(is),1:2),chanSorts.(myPols{sigPols(is)}).subScores(chanSorts.(myPols{sigPols(is)}).idx~=allSigs(is),1:2));
        qualVect(3) = klGetFNv2(chanSorts.(myPols{sigPols(is)}).subScores(chanSorts.(myPols{sigPols(is)}).idx==allSigs(is),1:2),chanSorts.(myPols{sigPols(is)}).subScores(chanSorts.(myPols{sigPols(is)}).idx~=allSigs(is),1:2),'sim',simMat);
        qualVect(4) = klGetFPv2(chanSorts.(myPols{sigPols(is)}).subScores(chanSorts.(myPols{sigPols(is)}).idx==allSigs(is),1:2),chanSorts.(myPols{sigPols(is)}).subScores(chanSorts.(myPols{sigPols(is)}).idx~=allSigs(is),1:2),'sim',simMat);
        spikes.qualVect = qualVect;
        
        if ~exist(sprintf('%s/%s/Channel%d/Unit%d/',myDir,file,ic,goodUnits),'file'),
            mkdir(sprintf('%s/%s/Channel%d/Unit%d/',myDir,file,ic,goodUnits));
        end
        save(sprintf('%s/%s/Channel%d/Unit%d/Spikes.mat',myDir,file,ic,goodUnits),'spikes');
        for ib = 1:length(sprintf('Getting Spiking Statistics...')), fprintf('\b'); end

    end
end
    fprintf('Channel %d completed\n',ic);
end
