function [epochMat, epochName, epochSlopes] = klParseEpochs(visSpks,visTimes,movSpks,movTimes,varargin)

% Set defaults
doSlope = 0;

preVis = [-100:0];
visTrans = [50:100];
visSust = [100:150];
preMov = [-50:0];
postMov = [0:50];
nextVis = [50:100];

epochs = {preVis,visTrans,visSust,preMov,postMov,nextVis};
spks  = {visSpks,visSpks,visSpks,movSpks,movSpks,movSpks};
times = {visTimes,visTimes,visTimes,movTimes,movTimes,movTimes};

epochName = {'Pre-Stim','VisTransient','VisSustained','Pre-Sacc','Post-Sacc','NextTrans'};

epochMat = cell2mat(cellfun(@getSpksFun,spks,times,epochs,'UniformOutput',0));

% for ir = 1:size(visSpks,1),
%     x = cellfun(@
%     epochMat(ir,1) = 

function [outMat, outSlope]=getSpksFun(spks,times,epochs,doSlope)
    out=nanmean(spks(:,ismember(times,epochs)),2);
    
    if doSlope,
        keyboard
    end
end

end