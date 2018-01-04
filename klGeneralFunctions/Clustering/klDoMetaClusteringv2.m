function [sortIDs, idxDist, raw, respSumStruct,rawLink] = klDoMetaClusteringv2(goodSDF,allTimeCell,varargin)

%% Set options and defaults
% restrict = 0;
zResp = [1];
setK = 5;
minN = 15;
metaN = 10;
show = 0;
unifK = 0;
randReps = 5;
meanMedian = 'median';
distOnly = 1;

% Default analysis epochs
preVis = -100:0;
visTrans = 50:100;
visSust = 100:150;
preMov = -50:0;
postMov = 0:50;
nextVis = 50:100;
myEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis};
myEpocInds = [1,1,1,2,2,2];
myEpocWinds = {[-200:300],[-300:200]};
normSD = 1:length(goodSDF);

%% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-k','k'}
            setK = varargin{varStrInd(iv)+1};
        case {'-n'}
            minN = varargin{varStrInd(iv)+1};
        case {'-mn'}
            metaN = varargin{varStrInd(iv)+1};
        case {'-r'}
            randReps = varargin{varStrInd(iv)+1};
        case {'-in'}
            respSumStruct = varargin{varStrInd(iv)+1};
        case {'-u'}
            unifK = varargin{varStrInd(iv)+1};
        case {'-m'}
            meanMedian = varargin{varStrInd(iv)+1};
        case {'-e'}
            myEpocs = varargin{varStrInd(iv)+1};
        case {'-ei'}
            myEpocInds = varargin{varStrInd(iv)+1};
        case {'-er'}
            myEpocWinds = varargin{varStrInd(iv)+1};
        case {'sd'}
            normSD = varargin{varStrInd(iv)+1};
            
    end
end

% Here are the sets of options/combinations to use
respNorms = {'ztr','ztrbl','max','none','bl','zbl'};
respSims = {'corr','euc'};%,'prod'};
respRand = {'pca'};
respInclude = {'mean','combo','slope','whole'};
zResp = 0;
doWhole = 1;

whichNorms = respNorms;
whichSims = respSims;
whichRand = respRand;
whichInclude = respInclude;
whichZ = [0,1];
whichWhole = 0;


if ~exist('respSumStruct','var'),
    %% Initialize output and printing variables
    isSameGroup = zeros(size(goodSDF{1},1),size(goodSDF{1},1),length(respNorms),length(respSims),length(respRand),length(respInclude),length(zResp));
    nNeeded = prod([length(respNorms),length(respSims),length(respRand),length(respInclude),length(zResp)]);
    nDone = 0;
    clear respSumStruct
    respSumStruct(1:nNeeded) = struct('ids',nan,'gap',nan,'gapErr',nan,'normType',nan,'simType',nan,'randType',nan,'include',nan,'k',nan);

    %% Loop on through and try the clustering procedures
    for in = 1:length(respNorms),
        for is = 1:length(respSims),
            for ir = 1:length(respRand),
                for ii = 1:length(respInclude),
                    for iz = 1:length(zResp),
                        clc;
                        nDone = nDone+1;
                        fprintf('Clustering combination %d of %d...\n',nDone,nNeeded);
                        [respK,respClustIDs,respGap,respGapErr,linkMat,linkInd, distMat] = klRespClustForMetaclustv3(goodSDF,allTimeCell,'norm',respNorms{in},'sim',respSims{is},'rand',respRand{ir},'resp',respInclude{ii},'z',zResp(iz),'-n',minN,'-r',randReps,'-e',myEpocs,'ei',myEpocInds,'rw',myEpocWinds,'sd',normSD,'do',distOnly);
                        respSumStruct(nDone).ids        = respClustIDs;
                        respSumStruct(nDone).gap        = respGap;
                        respSumStruct(nDone).gapErr     = respGapErr;
                        respSumStruct(nDone).normType   = respNorms{in};
                        respSumStruct(nDone).simType    = respSims{is};
                        respSumStruct(nDone).randType   = respRand{ir};
                        respSumStruct(nDone).include    = respInclude{ii};
                        respSumStruct(nDone).k          = respK;%size(respClustIDs,2);%setK;
                        respSumStruct(nDone).z          = zResp(iz);
                        respSumStruct(nDone).doWhole    = 0;
                        respSumStruct(nDone).linkMat    = linkMat;
                        respSumStruct(nDone).linkInd    = linkInd;
                        respSumStruct(nDone).distMat    = distMat;
                        if size(respClustIDs,2) >= setK,
                            for iii = 1:size(respClustIDs,1),
                                isSameGroup(iii,respClustIDs(:,setK)==respClustIDs(iii,setK),in,is,ir,ii,iz) = 1;
                            end
                        else
                            isSameGroup(:,:,in,is,ir,ii,iz) = nan;
                        end
                    end
                end
            end
        end
    end
end

%% Pull out the relevant procedures
includes = {respSumStruct.include};
norms = {respSumStruct.normType};
rands = {respSumStruct.randType};
sims = {respSumStruct.simType};
zs  = [respSumStruct.z];
wholes = [respSumStruct.doWhole];
includeParams = find(ismember(includes,whichInclude) & ismember(norms,whichNorms) & ismember(rands,whichRand) & ismember(sims,whichSims) & ismember(zs,whichZ) & ismember(wholes,whichWhole));

% Count up how often the clusters overlap
sumSame = zeros(size(goodSDF{1},1),size(goodSDF{1},1));
% sumSame = zeros(size(goodSDF,1),size(goodSDF,1));
allDist = nan(size(goodSDF{1},1),size(goodSDF{1},1),length(includeParams));
for ii = 1:length(includeParams)
    allDist(:,:,ii) = (respSumStruct(includeParams(ii)).distMat-nanmean(respSumStruct(includeParams(ii)).distMat(:)))./nanstd(respSumStruct(includeParams(ii)).distMat(:));
end

% Cluster the new distance matrix
switch meanMedian
    case {'mean'}
        raw = nanmean(allDist,3);
    case {'median'}
        raw = nanmedian(allDist,3);
end
[sortIDs,idxDist,~,~,rawLink] = klDistMatAgglomv2(raw,1:30,'-n',metaN);
while sum(isnan(sortIDs(:,end)))==size(sortIDs,1)
    sortIDs = sortIDs(:,1:(end-1));
end
