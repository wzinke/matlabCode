function [idx, outMap, outCount, outInds, outSim] = klAgglom_newVectv5(simMat,varargin)

% Set defaults
zeroTol = .0001;
cutLoopsN = 1000;
k = 1:5;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-k'},
            k = varargin{varStrInd(iv)+1};
    end
end

% Make vectors
simTri = triu(simMat);
simTri(simTri < zeroTol) = nan;
[sortSim,sortInds] = sort(simTri(:)); clear simTri;

colMat = repmat(1:size(simMat,2),size(simMat,1),1);
% rowMat = repmat((1:size(simMat,1))',1,size(simMat,2));
% sortRows = mod(((1:(size(simMat,1))^2)'-1),size(simMat,1))+1;
% sortRows = sortRows(sortInds(isfinite(sortSim)));
rowMat = colMat';
sortCols = colMat(sortInds(isfinite(sortSim))); clear colMat;
sortRows = rowMat(sortInds(isfinite(sortSim))); clear rowMat;
sortSim = sortSim(isfinite(sortSim));

% Initialize Outputs
outMap = uint16(zeros(size(simMat)));
outMap(:,1) = 1:length(simMat);
outSim = nan(1,length(simMat));

% Start loop
nLoopsAll = 0;
nLoopsGood = 0;
nLoopsSub = 1;
while nLoopsGood < (size(simMat,1)-1),%length(sortSim) > 1,

    % Advance total loop count (sortRows/sortCols/sortSim index)
    nLoopsAll = nLoopsAll+ 1;
    
    % Get the rows corresponding to the closest clusters
    rClust = outMap(sortRows(nLoopsAll),nLoopsGood+1);
    cClust = outMap(sortCols(nLoopsAll),nLoopsGood+1);
    
    if rClust == cClust,
        nLoopsSub = nLoopsSub+1;
        if nLoopsSub > cutLoopsN,
%             keyboard
            sortRows = outMap(sortRows(nLoopsAll:end),nLoopsGood+1);
            sortCols = outMap(sortCols(nLoopsAll:end),nLoopsGood+1);
            sortSim = sortSim(nLoopsAll:end);
            cutInds = sortRows==sortCols;
            if sum(cutInds) == length(sortSim),
                keyboard
            end
            sortRows(cutInds) = [];
            sortCols(cutInds) = [];
            sortSim(cutInds) = [];
            nLoopsAll = 0;
        end
        continue;
    else
        nLoopsGood = nLoopsGood + 1;
        nLoopsSub = 1;
    end
   
    % Save this similarity measurement
    outSim(nLoopsGood) = sortSim(nLoopsAll);
    
    % Get minimum (new cluster ID)
    newClust = min([rClust,cClust]);
    oldClust = max([rClust,cClust]);
    
    % Assign new members
    outMap(:,nLoopsGood+1) = outMap(:,nLoopsGood);
    outMap(outMap(:,nLoopsGood)==oldClust,nLoopsGood+1) = newClust;
end

% Now make an easy to access summary/companion matrices
% outCount = uint16(nan(size(outMap)));
% outInds = uint16(nan(size(outMap)));
% for i = 1:length(outMap),
%     uClusts = single(unique(outMap(:,i)));
%     [outCount(1:length(uClusts),i), sortInds] = sort(hist(single(outMap(:,i)),uClusts),'descend');
%     outInds(1:length(uClusts),i) = uClusts(sortInds);
% end
outCount = hist(outMap,1:length(outMap));
[outCount,outInds] = sort(outCount,'descend');

idx = uint8(ones(size(simMat,1),length(k)).*31);
for ik = 1:length(k),
    idxCol(ik) = find(outCount(k(ik),:)==max(outCount(k(ik),:)),1,'last');
    colKs = outInds(1:k(ik),idxCol(ik));
    for iik = 1:length(colKs),
        idx(outMap(:,idxCol(ik))==colKs(iik),ik) = iik;
    end
end

% keyboard