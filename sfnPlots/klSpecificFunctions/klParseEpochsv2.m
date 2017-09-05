function [epochMat, epochName, epochSlopes] = klParseEpochsv2(visSpks,visTimes,movSpks,movTimes,varargin)

% Set defaults
doSlope = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-s','slope'},
            doSlope = varargin{varStrInd(iv)+1};
    end
end

% Set time bins
preVis = [-100:0];
visTrans = [50:100];
visSust = [100:150];
preMov = [-50:0];
postMov = [0:50];
nextVis = [50:100];

epochs = {preVis,visTrans,visSust,preMov,postMov,nextVis};
spks  = {visSpks,visSpks,visSpks,movSpks,movSpks,movSpks};
times = {visTimes,visTimes,visTimes,movTimes,movTimes,movTimes};
slopes = {doSlope,doSlope,doSlope,doSlope,doSlope,doSlope};

epochName = {'Pre-Stim','VisTransient','VisSustained','Pre-Sacc','Post-Sacc','NextTrans'};

epochMat    = cell2mat(cellfun(@getSpksFun,spks,times,epochs,'UniformOutput',0));
epochSlopes = cell2mat(cellfun(@getSlopesFun,spks,times,epochs,'UniformOutput',0));

% for ir = 1:size(visSpks,1),
%     x = cellfun(@
%     epochMat(ir,1) = 

    function [outMat, outSlope]=getSpksFun(spks,times,epochs,doSlope)
        outMat=nanmean(spks(:,ismember(times,epochs)),2);
    end

    function outSlope = getSlopesFun(spks,times,epochs)
        [~,outSlope,~] = regression(repmat(epochs,size(spks,1),1),spks(:,ismember(times,epochs)));
    end
end