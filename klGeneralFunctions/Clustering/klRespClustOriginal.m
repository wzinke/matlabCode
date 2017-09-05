function [respK, respClustIDs, respGap, respGapErr, normVis, normMov] = klRespClustOriginal(allSDF,allTimes)

% global sdfK refine printProg respSimType nRespMembers refType kTypeResp nPrint respClustType pcaExplCrit pcaTypeWv pcaAgglomWv pcaWaveGapType pcaWaveRandType gapRandReps minGapConsec gapPerc gapDim
load('origClustParams.mat');

% respClustType = 'pca';
gapDim = 1;
doWeights = 1;
refine = 1;
pcaRefType = 'same';

if ~any(size(allSDF)==1),
    fprintf('\tAligning SDF and times...');
    visZero = cellfun(@(times) find(times == 0,1),allTimes(:,1));
    movZero = cellfun(@(times) find(times == 0,1),allTimes(:,2));
    rawZero = cellfun(@(times) find(times == 0,1),allTimes(:,3));
    visAlign = klAlignv2(allSDF(:,1),visZero);
    visTimes = nanmean(klAlignv2(allTimes(:,1),visZero),1);
    movAlign = klAlignv2(allSDF(:,2),movZero);
    movTimes = nanmean(klAlignv2(allTimes(:,2),movZero),1);
    rawAlign = klAlignv2(allSDF(:,3),visZero);
    rawTimes = nanmean(klAlignv2(allTimes(:,3),rawZero),1);

    visAlign(visAlign == inf) = nan;
    movAlign(movAlign == inf) = nan;
    rawAlign(rawAlign == inf) = nan;
    fprintf('Done!\n');
else
    visAlign = allSDF{1}; visTimes = allTimes{1};
    movAlign = allSDF{2}; movTimes = allTimes{2};
    rawAlign = allSDF{3}; rawTimes = allTimes{3};
end

nRespMembers = 15;
respSimType = 'corr';
respValType = 'means';
respNormType = 'z';

switch respNormType,
    case {'z'},
        blMeans = nanmean(visAlign(:,ismember(visTimes,-300:-200)),2);
        blStds = nanstd(visAlign(:,ismember(visTimes,-300:-200)),[],2);
        normVis = (visAlign-repmat(blMeans,1,size(visAlign,2)))./repmat(blStds,1,size(visAlign,2));
        normMov = (movAlign-repmat(blMeans,1,size(movAlign,2)))./repmat(blStds,1,size(movAlign,2));
    case {'max'};
        blMeans = nanmean(visAlign(:,ismember(visTimes,-300:-200)),2);
        catVals = [visAlign(:,ismember(visTimes,-200:300)),movAlign(:,ismember(movTimes,-300:200))];
        absMax = max(abs(catVals-repmat(blMeans,1,size(catVals,2))),[],2);
        normVis = (visAlign-repmat(blMeans,1,size(visAlign,2)))./repmat(absMax,1,size(visAlign,2));
        normMov = (movAlign-repmat(blMeans,1,size(movAlign,2)))./repmat(absMax,1,size(movAlign,2));
    case {'none'},
        normVis = visAlign;
        normMov = movAlign;
end
        
catNorms = [normVis(:,ismember(visTimes,-200:300)),normMov(:,ismember(movTimes,-300:200))];
goodRows = sum(isfinite(catNorms),2) == size(catNorms,2);

preVis = [-100:0];
visTrans = [50:100];
visSust = [100:150];
preMov = [-50:0];
postMov = [0:50];
nextVis = [50:100];
allEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis};

normVis = normVis(goodRows,:);
normMov = normMov(goodRows,:);

[respMat, respNames, respSlope] = klParseEpochsv2(normVis, visTimes, normMov, movTimes,'-s',1);
% respMat = [respMat,respSlope];

switch respValType,
    case {'means','vals'},
        respIn = respMat;
    case {'slope','slopes','s'},
        respIn = respSlope;
    case {'combined'},
        respIn = [respMat,respSlope];
    otherwise
        respIn = respMat;
end
respIn = (respIn-repmat(nanmean(respIn,1),size(respIn,1),1))./repmat(nanstd(respIn,[],1),size(respIn,1),1);

switch respClustType
    case 'agglom'
%         [visClustIDs, visClusts, visSimMat, visLink] = klAgglomv10(visAlign(:,ismember(nanmean(visTimes,1),vWind)),'-k',sdfK,'-r',refine,'-p',printProg,'-t',visSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeVis,'nprint',nPrint);
%         [movClustIDs, movClusts, movSimMat, movLink] = klAgglomv10(movAlign(:,ismember(nanmean(movTimes,1),mWind)),'-k',sdfK,'-r',refine,'-p',printProg,'-t',movSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeMov,'nprint',nPrint);
        [respClustIDs, respClusts, respSimMat, respLink] = klAgglomv10(respIn,'-k',sdfK,'-r',refine,'-p',printProg,'-t',respSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeResp,'nprint',nPrint);
%         while sum(isnan(visClustIDs(:,end))) == size(visClustIDs,1),
%             visClustIDs = visClustIDs(:,1:(end-1));
%         end
%         while sum(isnan(movClustIDs(:,end))) == size(movClustIDs,1),
%             movClustIDs = movClustIDs(:,1:(end-1));
%         end
        while sum(isnan(respClustIDs(:,end))) == size(respClustIDs,1),
            respClustIDs = respClustIDs(:,1:(end-1));
        end
    case 'kmeans'
        warning off
        tmpVis = visAlign(:,ismember(nanmean(visTimes,1),visWind));
        tmpMov = movAlign(:,ismember(nanmean(movTimes,1),movWind));
        for ir = 1:size(visAlign,1),
            firstValV(ir) = find(~isnan(tmpVis(ir,:)),1);
            lastValV(ir) = find(~isnan(tmpVis(ir,:)),1,'last');
            firstValM(ir) = find(~isnan(tmpMov(ir,:)),1);
            lastValM(ir) = find(~isnan(tmpMov(ir,:)),1,'last');
        end
        clpVis = tmpVis(:,max(firstValV):min(lastValV));
        clpMov = tmpMov(:,max(firstValM):min(lastValM));
        for ik = 1:length(sdfK),
            for ir = 1:nReps,
                [icdV(:,ir)] = kmeans(clpVis,sdfK(ik),'Distance',kTypeResp);
                [icdM(:,ir)] = kmeans(clpMov,sdfK(ik),'Distance',kTypeResp);
            end
            visClustIDs(:,ik) = klModeKmeans(icdV);
            movClustIDs(:,ik) = klModeKmeans(icdM);
            clear icdV icdM;
        end
        visSimMat = [];
        movSimMat = [];
    case 'pca'
        warning off
        [~,scores,~,~,expl] = pca(respIn);
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
                [respClustIDs, respClusts, respSimMat, respLink] = klAgglomv10(scores,'-k',sdfK,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nRespMembers,'refType',pcaRefType,'kType',kTypeResp,'nprint',nPrint,'-w',weights);
            otherwise
                [respClustIDs, respClusts, respSimMat, respLink] = klAgglomv10(scores,'-k',sdfK,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nRespMembers,'refType',pcaRefType,'kType',kTypeResp,'nprint',nPrint,'-w',weights);
        end
        while sum(isnan(respClustIDs(:,end))) == size(respClustIDs,1),
            respClustIDs = respClustIDs(:,1:(end-1));
        end
    otherwise
        [respClustIDs, respClusts, respSimMat, respLink] = klAgglomv10(respIn,'-k',sdfK,'-r',refine,'-p',printProg,'-t',respSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeResp,'nprint',nPrint);
        while sum(isnan(respClustIDs(:,end))) == size(respClustIDs,1),
            respClustIDs = respClustIDs(:,1:(end-1));
        end
end
% [visSimF,visSimR] = klGetF(visAlign(:,ismember(nanmean(visTimes,1),visWind)),visClustIDs,'-t',visSimType,'-s',visSimMat);
% [movSimF,movSimR] = klGetF(movAlign(:,ismember(nanmean(movTimes,1),movWind)),movClustIDs,'-t',movSimType,'-s',movSimMat);
% [respSimF, respSimR] = klGetF(respMat,respClustIDs,'-t',respSimType,'-s',respSimMat);

% [visK,visGap,visGapErr] = klGetGapv10(visAlign(:,ismember(nanmean(visTimes,1),visWind)),visClustIDs,'-t',visSimType,'-r',visGapType,'randval',visRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim);
% [movK,movGap,movGapErr] = klGetGapv10(movAlign(:,ismember(nanmean(movTimes,1),movWind)),movClustIDs,'-t',movSimType,'-r',movGapType,'randval',movRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim);
[respK,respGap,respGapErr] = klGetGapv10(respIn,respClustIDs,'-t',respSimType,'-r',respGapType,'randval',respRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim,'refine',refine);

respClustTmp = respClustIDs;
normVisTmp = normVis; normMovTmp = normMov;
clear respClustIDs normVis normMov;
respClustIDs = nan(size(catNorms,1),size(respClustTmp,2));
normVis = nan(size(catNorms,1),size(normVisTmp,2));
normMov = nan(size(catNorms,1),size(normMovTmp,2));
respClustIDs(goodRows,:) = respClustTmp;
normVis(goodRows,:) = normVisTmp;
normMov(goodRows,:) = normMovTmp;
end