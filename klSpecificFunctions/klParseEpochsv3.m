function [epochMat, epochSlopes] = klParseEpochsv3(spks,times,epochs,varargin)

epochMat    = cell2mat(cellfun(@getSpksFun,spks,times,epochs,'UniformOutput',0));
epochSlopes = cell2mat(cellfun(@getSlopesFun,spks,times,epochs,'UniformOutput',0));

    function outMat=getSpksFun(spks,times,epochs)
        outMat=nanmean(spks(:,ismember(times,epochs)),2);
    end

    function outSlope = getSlopesFun(spks,times,epochs)
        [~,outSlope,~] = regression(repmat(epochs,size(spks,1),1),spks(:,ismember(times,epochs)));
    end
end