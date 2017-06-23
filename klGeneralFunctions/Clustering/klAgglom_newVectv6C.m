function [idx, outMap, outCount, outInds, outSim] = klAgglom_newVectv6C(simMat,varargin)

% Set defaults
zeroTol = .0001;
cutLoopsN = 1000;
k = 1:5;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-k'}
            k = varargin{varStrInd(iv)+1};
    end
end

% Make vectors
simTri = triu(simMat);
simTri(simTri < zeroTol) = nan;
% No need for simMat except for its dims
sizeSimMat=size(simMat); clear simMat;
[sortSim,sortInds] = sort(simTri(:)); clear simTri;
% Prune nans and index of nans
firstNan=find(isnan(sortSim),1,'first');
sortSim=sortSim(1:firstNan-1);
sortInds=sortInds(1:firstNan-1);
% Convert ind to r,c
[sortRows,sortCols]=ind2sub(sizeSimMat,sortInds);

% Initialize Outputs
outMap = uint16(zeros(sizeSimMat));
outMap(:,1) = 1:sizeSimMat(1);
outSim = nan(1,sizeSimMat(1));

% Start loop
nLoopsAll = 0;
nLoopsGood = 0;
nLoopsSub = 1;
while nLoopsGood < sizeSimMat(1)-1 %length(sortSim) > 1,
    
    % Advance total loop count (sortRows/sortCols/sortSim index)
    nLoopsAll = nLoopsAll+ 1;
    
    % Get the rows corresponding to the closest clusters
    rClust = outMap(sortRows(nLoopsAll),nLoopsGood+1);
    cClust = outMap(sortCols(nLoopsAll),nLoopsGood+1);
    
    if rClust == cClust
        nLoopsSub = nLoopsSub+1;
        if nLoopsSub > cutLoopsN
            %             keyboard
            sortRows = outMap(sortRows(nLoopsAll:end),nLoopsGood+1);
            sortCols = outMap(sortCols(nLoopsAll:end),nLoopsGood+1);
            sortSim = sortSim(nLoopsAll:end);
            cutInds = sortRows==sortCols;
            if sum(cutInds) == length(sortSim)
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

% Comment: At this point Vars needed are
% output -> idx, outMap, outCount, outInds, outSim
% vars are:
%   cClust                1x1                      2  uint16
%   cutInds         2259440x1                2259440  logical
%   cutLoopsN             1x1                      8  double
%   iv                    1x1                      8  double
%   k                     1x1                      8  double
%   nLoopsAll             1x1                      8  double
%   nLoopsGood            1x1                      8  double
%   nLoopsSub             1x1                      8  double
%   newClust              1x1                      2  uint16
%   oldClust              1x1                      2  uint16
%   outCount           3000x3000            72000000  double
%   outInds            3000x3000            72000000  double
%   outMap             3000x3000            18000000  uint16
%   outSim                1x3000               24000  double
%   rClust                1x1                      2  uint16
%   simMat             3000x3000            72000000  double
%   sortCols        2250000x1                4500000  uint16
%   sortInds        9000000x1               72000000  double
%   sortRows        2250000x1                4500000  uint16
%   sortSim         2250000x1               18000000  double
%   varStrInd             1x1                      8  double
%   varargin              1x2                    236  cell
%   zeroTol               1x1                      8  double
%
% we only need the following vars from here on...>
%   outCount           3000x3000            72000000  double
%   outInds            3000x3000            72000000  double
%   outMap             3000x3000            18000000  uint16
%   outSim                1x3000               24000  double
% need sizeSimMat_1=size(simMat,1)
% k
t=whos;
disp(['Bytes before cleaning: ',num2str(sum([t.bytes]))]);
% so that we can clear all vars not begin with out...
outSizeSimMat_1=sizeSimMat(1);
outK=k;
% clean not needed
clear -regexp '^((?!out).)*$'
t=whos;
disp(['Bytes **After** cleaning: ',num2str(sum([t.bytes]))]);


% Now make an easy to access summary/companion matrices
% very nice ...
outCount = hist(outMap,1:length(outMap));
[outCount,outInds] = sort(outCount,'descend');

idx = uint8(ones(outSizeSimMat_1,length( outK)).*31);
for ik = 1:length(outK)
    idxCol(ik) = find(outCount(outK(ik),:)==max(outCount(outK(ik),:)),1,'last');
    colKs = outInds(1:outK(ik),idxCol(ik));
    for iik = 1:length(colKs)
        idx(outMap(:,idxCol(ik))==colKs(iik),ik) = iik;
    end
end

% keyboard