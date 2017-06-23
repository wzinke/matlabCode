function [idx, outMap,outCount,outInds, outSim] = klAgglom_newv5a(simAll,varargin)

% Set defaults
linkType = 'average';
nPrint = 1000;
print = 1;
k = 1:6;
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)},
        case {'-p','print'},
            print = varargin{varStrInd(iv)+1};
        case {'-k','k'},
            k = varargin{varStrInd(iv)+1};
    end
end

% Save raw input matrix
totalN = size(simAll,1);
if isnan(simAll(2,1)),
    for ir = 1:size(simAll,1),
        for ic = ir:size(simAll,1),
            if ir ~= ic,
                simAll(ic,ir) = simAll(ir,ic);
            end
        end
    end
end
if ~isnan(simAll(1,1)),
    for ir = 1:size(simAll,1),
        simAll(ir,ir) = nan;
    end
end

% Initialize Output
clustMembs = cell(size(simAll,1),size(simAll,2));
clustMembs(:,1) = mat2cell((1:size(simAll,1))',ones(size(simAll,1),1));
outSim = nan(1,size(simAll,1));
clustMap = 1:size(simAll,1);

% Start the loop
maxWhiles = 0;
simCutT = nan(1,size(simAll,1));
assignT = nan(1,size(simAll,1));
outMap = nan(size(simAll,1));
outMap(:,1) = 1:size(simAll,1);

% Keep a counter of which original rows remain valid
goodCols = 1:size(simAll,1);
for ii = 1:(size(simAll,1)-1),
    fulltic = tic;
    if print && mod(ii,nPrint)==0,
        fprintf('Doing loop %d of %d\n',ii,totalN);
    end
    % Get the minimum similarity
    findtic = tic;
    [minRow,minCol] = find(simAll==min(simAll(:)),1);
    % Take out the minimum entry from simMat
    outSim(ii) = simAll(minRow,minCol);
    nLoops =0;
    findTime(ii) = toc(findtic);
    
    % Get current cluster members
    % clustMap points to where that row's current identity lies
%     rowMembs = clustMembs{minRow,ii};
%     colMembs = clustMembs{minCol,ii};
%     rowMembs = clustMembs{clustMap(minRow),ii};
%     colMembs = clustMembs{clustMap(minCol),ii};
    rowMembs = find(outMap(:,ii)==outMap(minRow,ii));
    colMembs = find(outMap(:,ii)==outMap(minCol,ii));
    
    % Get new cluster members
    newMembs = unique([rowMembs;colMembs]);
    
    x=tic;
    minInd = min(clustMap([minRow,minCol]));
    maxInd = max(clustMap([minRow,minCol]));
    
    % Change up simMat and clustMembs
    switch linkType
        case {'single'},
            % Place minimum distances for this group into the first
            % position
            newVals = min([simAll(minRow,:);simAll(:,minCol)'],[],1);
%             newVals(max([minRow,minCol])) = [];
%             newVals(max([minRow,minCol])) = nan;
            newVals(newMembs) = nan;
        case {'average'}
            newVals = nansum([(simAll(minRow,:).*length(rowMembs));(simAll(:,minCol)'.*length(colMembs))],1)./length(newMembs);
%             newVals(max([minRow,minCol])) = nan;
            newVals(newMembs) = nan;
        case {'sum'},
            newVals = nansum([simAll(minRow,:);simAll(:,minCol)'],1);
            newVals(max([minRow,minCol])) = [];
        case {'centroid'}
    end
    newValTime(ii) = toc(x);
    
    cutMatTic = tic;
%     simAll = simAll(~ismember(1:size(simAll,1),max([minRow,minCol])),~ismember(1:size(simAll,1),max([minRow,minCol])));
    tic1 = tic;
%     simAll = simAll([1:(maxInd-1),(maxInd+1):end],[1:(maxInd-1),(maxInd+1):end]);
    simAll(newMembs,newMembs) = deal(nan);
    simAll(minInd,:) = newVals;
    simAll(:,minInd) = newVals;
    cutMatTime(ii) = toc(cutMatTic);
    
    outMap(:,ii+1) = outMap(:,ii);
    outMap(newMembs,ii+1) = minInd;
    
    % Put newMembs in both old spots so they can be accessed either
    % way
    y=tic;
%     clustMembs(goodCols,ii+1) = clustMembs(goodCols,ii);
    clustMembs(:,ii+1) = clustMembs(:,ii);
    clustMembs{minInd,ii+1} = newMembs;
    clustMembs{maxInd,ii+1} = [];
%     [clustMembs{newMembs(2:end),ii+1}] = deal([]);
    % Update clustMap -> point to the row on line above when other
    % newMembs members are needed (this should clear up memory and
    % help with finding clusters after the fact
    assignT(ii) = toc(y);
    clustMap(newMembs) = newMembs(1);
    fullTime(ii) = toc(fulltic);
end
    
% Now make an easy to access summary/companion matrices
outCount = uint16(nan(size(outMap)));
outInds = uint16(nan(size(outMap)));
for i = 1:length(outMap),
    uClusts = single(unique(outMap(:,i)));
    [outCount(1:length(uClusts),i), sortInds] = sort(hist(single(outMap(:,i)),uClusts),'descend');
    outInds(1:length(uClusts),i) = uClusts(sortInds);
end

%     
%     keyboard
    
idx = uint8(ones(size(simAll,1),length(k)).*31);
for ik = 1:length(k),
    idxCol(ik) = find(outCount(k(ik),:)==max(outCount(k(ik),:)),1,'last');
    colKs = outInds(1:k(ik),idxCol(ik));
    for iik = 1:length(colKs),
        idx(outMap(:,idxCol(ik))==colKs(iik),ik) = iik;
    end
end

    
    
    
    
    
    
    