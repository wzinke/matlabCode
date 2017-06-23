function [outK, idx, idxCol, aggMap, aggCount, aggInds, aggSim, randStruct] = klAgglomClust_newv4(obs,varargin)

% Set defaults
distType = 'euc';
randReps = 10;
distros = ones(1,size(obs,2));
minMembs = .01;
sizeDir =  'last';
k = 1:6;
qualType = 'cross';
matchKs = 0;
tempFold = '.';
saveFlag = '';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)},
        case {'-d','dist'},
            distType = varargin{varStrInd(iv)+1};
        case {'-distro'},
            distros = varargin{varStrInd(iv)+1};
        case {'-m','min'},
            minMembs = varargin{varStrInd(iv)+1};
	case {'-s'},
	    saveFlag = varargin{varStrInd(iv)+1};
	case {'-t'},
	    tempFold = varargin{varStrInd(iv)+1};
    end
end

distMat = klGetDist_newv1(obs,distType);
distMatGonly = klGetDist_newv1(obs(:,distros==1),distType);

% Get random parameters
nObs = size(obs,1); 
mus = nanmean(obs,1);
sigs = nanstd(obs,[],1);
mins = min(obs,[],1);
maxs = max(obs,[],1);
par1(distros==1) = mus(distros==1);
par2(distros==1) = sigs(distros==1);
par1(distros==2) = mins(distros==2);
par2(distros==2) = maxs(distros==2);
clear obs;

% Cluster the actual observations
aggTic = tic;
[idx, aggMap, aggCount, aggInds, aggSim] = klAgglom_newVectv5(distMat,'-k',k);

%[aggMap, aggCount, aggInds] = klAgglom_newv5a(distMat);
aggTime = toc(aggTic);
fprintf('Actual observation finished in %s\n',printTiming(aggTic));

% Save and clear for space
save(sprintf('%s/tempClustSave%s.mat',tempFold,saveFlag),'aggMap','aggCount','aggInds','aggSim','-v7.3');
aggMap = []; aggCount = []; aggInds = []; aggSim = [];

% Loop through K to get a cluster quality measure...
% Maybe maximum intra-cluster distance?
if minMembs < 1,
    minMembs = floor(minMembs.*nObs);
end

% Get random observaitons

% Randomize and "sort"
randMax = nan(nObs,randReps);
fullRandTic = tic;
randT = nan(1,randReps);
for ir = 1:randReps,
    randTic = tic;
    % Make a random observation and get distance matrices
    randObs = klMakeRand(nObs,par1,par2,'-d',distros);
    randMat = klGetDist_newv1(randObs,distType);

    % Do the clustering
    [~, randCount, ~] = klAgglom_newVectv5(randMat);
    %[~, randCount, ~] = klAgglom_newv5a(randMat);
    randMax(:,ir) = max(randCount,[],2);
    clear randCount randObs randMat
    randT(ir) = toc(randTic);
    fprintf('Random observation %d finished in %s\n',ir, printTiming(randTic));
end
fullRandT = toc(fullRandTic);
fprintf('All random observations finished in %s\n',printTiming(fullRandTic));

% [randMap, randCount, randInds] = klRandObsPar(nObs,par1,par2,distros,randReps,distType);

% % Get mean and standard deviation of the random observation counts
% blankCell = cell(1,randReps);
% dimCell = cell(1,randReps); [dimCell{:}] = deal(2);
% maxes = cellfun(@max,randCount,blankCell,dimCell);

maxRandMems = nanmean(randMax(k,:),2)+2*nanstd(randMax(k,:),[],2);
randStruct.max = randMax;
randStruct.maxMembs = maxRandMems;
randStruct.par1 = par1;
randStruct.par2 = par2;
randStruct.distros = distros;

clearvars -except randStruct k matchKs maxRandMembs nObs tempFold saveFlag idx

hasK = nan(1,length(k));
colInd = nan(1,length(k));
load(sprintf('%s/tempClustSave%s.mat',tempFold,saveFlag));
if matchKs,
    for ik = 1:length(k),
        if any(aggCount(k(ik),:) > randStruct.maxMembs(k(ik))),
            hasK(k(ik)) = 1;
            colInd(k(ik)) = find(aggCount(k(ik),:) > randStruct.maxMembs(k(ik)),1,'last');
        else
            hasK(k(ik)) = 0;
            colInd(k(ik)) = nan;
        end
    end
else
    for ik = 1:length(k),
        if any(aggCount(k(ik),:) > randStruct.maxMembs(2)),
            hasK(k(ik)) = 1;
            colInd(k(ik)) = find(aggCount(k(ik),:) > randStruct.maxMembs(2),1,'last');
        else
            hasK(k(ik)) = 0;
            colInd(k(ik)) = nan;
        end
    end
end


outK = find(hasK,1,'last');

