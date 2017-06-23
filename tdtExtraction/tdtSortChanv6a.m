function [outSorts] = tdtSortChanv6a(chanData,varargin)

warning off

% Set defaults
maxWvs = 20000;
maxAgglom = 2500;
minRefract = .6;
% thisThresh = -4*(median(abs(chanData)/.6745));
thisThresh = 2*rms(chanData);
spkFreq = 1000/24414;
wvTimes = ((1:32)-9).*spkFreq;
chanTimes = single(0:spkFreq:(spkFreq*(size(chanData,2)-1))); % 1000x multiplier converts to ms
compType = 'euc';
doPos = 0;
smooth = 0;
align = 0;
nDims = 2;
dimRedType = 'pca';
doTime = 0;
distros = ones(1,nDims);
minSNR = 1.25;

% Decode varargin
varStrInd=find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t'},
            thisThresh = varargin{varStrInd(iv)+1};
        case {'-p'},
            doPos = varargin{varStrInd(iv)+1};
        case {'-c'},
            compType = varargin{varStrInd(iv)+1};
        case {'-d','nDims'},
            nDims = varargin{varStrInd(iv)+1};
    end
end

% Set threshold, get crossing times
% [spkTimesN,spkWavesN,spkThreshN] = klThreshCrossv5(chanData,'times',chanTimes,'-m',minRefract,'-t',-thisThresh);
% [spkTimesP,spkWavesP,spkThreshP] = klThreshCrossv5(chanData,'times',chanTimes,'-m',minRefract,'-t',thisThresh,'-p',1);

[spkTimes,spkWaves,spkPols] = klThreshCrossv6(chanData,'times',chanTimes,'-m',minRefract,'-t',thisThresh);
% clear chanData
chanData = [];
spkTimesN = spkTimes(spkPols==2);
spkTimesP = spkTimes(spkPols==1); spkTimes = [];
spkWavesN = spkWaves(spkPols==2,:);
spkWavesP = spkWaves(spkPols==1,:); spkWaves = [];


outSorts.neg = tdtSortWavesv1a(spkTimesN,spkWavesN,'-d',nDims);
outSorts.pos = tdtSortWavesv1a(spkTimesP,spkWavesP,'-d',nDims);