function [SNR,isoScore,fnScore,fpScore] = klGetAllIso(waves,sortStruct,varargin)

upSamp = 0;
times = 1:size(waves,2);
noiseC = 5;
nnPerc = .05;
faCut = .5;
sortType = 'kmeans';
distType = 'euc';
lamb = 10;
maxWaves = 10000;
silent = 1;
forceK = 0;
align = 1;
reSort = 1;
inK = 2;
wvTimes = ((1:32)-9).*25;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','times'},
            times = varargin{varStrInd(iv)+1};
        case {'-c','c'},
            noiseC = varargin{varStrInd(iv)+1};
        case {'-u','upsamp'},
            upSamp = varargin{varStrInd(iv)+1};
        case {'-a','align'},
            align = varargin{varStrInd(iv)+1};
        case {'-n','nnperc'}
            nnPerc = varargin{varStrInd(iv)+1};
        case {'type'},
            sortType = varargin{varStrInd(iv)+1};
        case {'-d'},
            distType = varargin{varStrInd(iv)+1};
        case {'-l'},
            lamb = varargin{varStrInd(iv)+1};
        case {'-k'},
            inK = varargin{varStrInd(iv)+1};
            forceK = 1;  
        case {'resort','-r'},
            if isStruct(varargin{varStrInd(iv)+1}),
                oldSortStruct = varargin{varStrInd(iv)+1};
                reSort = 0;
            else
                fprintf('Missing required input... resorting anyway\n');
            end
    end
end

% Get subset of waveforms
subWaves = waves(sortStruct.waveInds,:);

% Make them acceptable for spline interpolation
cutWaves = klCutForSplinesv1(subWaves);
splWaves = nan(size(cutWaves,1),length(1:.1:32));
for iw = 1:size(cutWaves,1),
    splWaves(iw,:) = spline(1:32,cutWaves(iw,:),1:.1:32);
end
[alWaves,alTimes] = klTroughAlignv5(splWaves,spline(1:32,wvTimes,1:.1:32),0);

% Get spike and noise clusters
spkClust = alWaves(sortStruct.idx == sortStruct.spkNum,:);
nzClust  = alWaves(sortStruct.idx ~= sortStruct.spkNum,:);

%% Start by getting SNR
SNR = klGetSNRv1(spkClust);

uID = unique(sortStruct.idx); uID(isnan(uID)) = [];
if length(uID) > 1,
    %% Get Isolation Score
    [isoScore, simMat] = klGetISv2(spkClust,nzClust);

    %% Let's work on false negative scores
    fnScore = klGetFNv2(spkClust,nzClust,'sim',simMat);

    %% Now false positive scores
    fpScore = klGetFPv2(spkClust,nzClust,'sim',simMat);
else
    isoScore = nan; fnScore = nan; fpScore = nan;
end
end
        