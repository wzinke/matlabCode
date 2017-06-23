function [wvK, waveClustIDs, gapVect, gapErr] = klWaveClustOriginal(waves)

% global k refine printProg waveSimType nWvMembers refType kTypeWv nPrint wvClustType pcaExplCrit pcaTypeWv pcaAgglomWv pcaWaveGapType pcaWaveRandType gapRandReps minGapConsec gapPerc gapDim
load('origClustParams.mat');

% nWvMembers
nWvMembers = 5;
doWeights = 1;
refine = 0;


% waveSimType = 'pca';

%% Cluster waveforms
startTic = tic;
fprintf('\n*** Clustering Waveforms ***\n');
switch wvClustType
    case 'agglom'
        [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(waves,'-k',k,'-r',refine,'-p',printProg,'-t',waveSimType,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint);
        wvSimInd = klGetF(waves,waveClustIDs,'-t',waveSimType,'-s',simMat);
    case 'kmeans'
        warning off
        for ir = 1:size(allWaves,1),
            firstVal(ir) = find(~isnan(waves(ir,:)),1);
            lastVal(ir) = find(~isnan(waves(ir,:)),1,'last');
        end
        clpWvs = waves(:,max(firstVal):min(lastVal));
        for ik = 1:length(k),
            for ir = 1:nReps,
                [icd(:,ir)] = kmeans(clpWvs,k(ik),'Distance',kTypeWv);
            end
            waveClustIDs(:,ik) = klModeKmeans(icd);
            clear icd;
        end
        simMat = [];
        wvSimInd = klGetF(waves,waveClustIDs,'-t',waveSimType,'-s',simMat);
    case 'pca'
        warning off
        for ir = 1:size(waves,1),
            firstVal(ir) = find(~isnan(waves(ir,:)),1);
            lastVal(ir) = find(~isnan(waves(ir,:)),1,'last');
        end
        clpWvs = waves(:,max(firstVal):min(lastVal));
        [~,scores,~,~,expl] = pca(clpWvs);
        cumExpl = cumsum(expl);
        if doWeights,
            weights = expl';
        else
            weights = ones(size(expl'));
        end
        nDims = find(cumExpl > pcaExplCrit,1);
        switch pcaTypeWv,
            case 'kmeans',
                for ik = 1:length(k),
                    for ir = 1:nReps,
                        [pcaIDX(:,ir),~,pcaDist(1:k(ik),ir)] = kmeans(scores(:,1:nDims),k(ik));
                    end
                    waveClustIDs(:,ik) = klModeKmeans(pcaIDX,pcaDist);
                end
                simMat = [];
            case 'agglom'
%                 [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(scores,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint,'-w',weights);
%                 [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(scores(:,1:nDims),'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint,'-w',weights(1:nDims));
                [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(scores,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint,'-w',weights);
            otherwise
                [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(scores,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',pcaRefType,'kType',kTypeWv,'nprint',nPrint,'-w',weights);
        end
        wvSimInd = klGetF(scores(:,1:nDims),waveClustIDs,'-t',pcaAgglomWv,'-s',simMat);
%                 [wvK, gapVect, gapErr] = klGetGapv10(waves,waveClustIDs,'-t',waveSimType);
    case 'lpp'
        warning off
        for ir = 1:size(allWaves,1),
            firstVal(ir) = find(~isnan(waves(ir,:)),1);
            lastVal(ir) = find(~isnan(waves(ir,:)),1,'last');
        end
        clpWvs = waves(:,max(firstVal):min(lastVal));
        [outVals,~,eigVals] = klLPPv1(clpWvs);
        if doWeights,
            weights = eigVals';
        else
            weights = ones(size(eigVals'));
        end
        
%                 cumExpl = cumsum(expl);
%                 nDims = find(cumExpl > pcaExplCrit,1);
        switch pcaTypeWv,
            case 'kmeans',
                for ik = 1:length(k),
                    for ir = 1:nReps,
                        [pcaIDX(:,ir),~,pcaDist(1:k(ik),ir)] = kmeans(outVals,k(ik));
                    end
                    waveClustIDs(:,ik) = klModeKmeans(pcaIDX,pcaDist);
                end
                simMat = [];
            case 'agglom'
                [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(outVals,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint,'-w',weights);
            otherwise
                [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(outVals,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',pcaRefType,'kType',kTypeWv,'nprint',nPrint,'-w',weights);
        end
        wvSimInd = klGetF(scores(:,1:nDims),waveClustIDs,'-t',pcaAgglomWv,'-s',simMat);
%                 [wvK, gapVect, gapErr] = klGetGapv10(waves,waveClustIDs,'-t',waveSimType);
    otherwise
        [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(waves,'-k',k,'-r',refine,'-p',printProg,'-t',waveSimType,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint);
        wvSimInd = klGetF(waves,waveClustIDs,'-t',waveSimType,'-s',simMat);
end
% Clear out columns of all-nan
while sum(isnan(waveClustIDs(:,end))) == size(waveClustIDs,1),
    waveClustIDs = waveClustIDs(:,1:(end-1));
end

fprintf('\n*** Waveforms Clustered in %s  ***\n',printTiming(startTic));

if strcmp(wvClustType,'pca'),
    [wvK, gapVect, gapErr] = klGetGapv12(scores,waveClustIDs,'-t',pcaAgglomWv,'-n',nWvMembers,'-r',pcaWaveGapType,'randval',pcaWaveRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim,'-w',weights,'refine',refine);
elseif strcmp(wvClustType,'lpp'),
    [wvK, gapVect, gapErr] = klGetGapv12(outVals,waveClustIDs,'-t',pcaAgglomWv,'-n',nWvMembers,'-r',pcaWaveGapType,'randval',pcaWaveRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim,'-w',weights,'refine',refine);
else
    [wvK, gapVect, gapErr] = klGetGapv12(waves,waveClustIDs,'-t',waveSimType,'-n',nWvMembers,'-r',wvGapType,'randval',wvRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim,'refine',refine);
end
 
