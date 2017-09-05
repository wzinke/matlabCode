% %  klUnitAgglom

%% Check current directory
if ~exist('./klUnitAgglomv13.m','file'),
    fprintf('Please change directory to root MATLAB folder...\n');
    keyboard
end

grandTic = tic;

close all;
reselect = 0;
if ~reselect, % && exist('allPseudoF','var')% && sum(isnan(allPseudoF(:))) ~= numel(allPseudoF),
    freshRun        = 0;
    reextract       = 1;
    reClust         = 1;
    reClustResp     = 0;
    reAlignSDF      = 0;
    if freshRun,
        clearvars -except grandTic;
        reClust     = 1;
        reextract   = 1;
        reClustResp = 1;
        reAlignSDF  = 1;
    end

    %% Pull in globals
    global masterSorts excelNum excelAll wvClusts catFacts factNames keepFacts rawWaves allWaves normWaves wvTimes monkID catCrit
    
    %% Set constants
    monk = {'Gauss','Helmholtz','Darwin'};
%     monk = {'Gauss'};
%     monk = {'Helmholtz'};
%     monk = {'Darwin'};
    
    % Subset selection
    task            = 'MG';
    areas           = {'FEF'};
    minWidth        = 0;
    cutWideNoise    = 1;
    isiPerc         = 1;
    isoMin          = 0;
    fpMax           = 1;
    snrMin          = 1;
    cutTails        = 0;
    sortedOnly      = 1;
    doCuts          = 1;
    
    % Clustering Constants
    k               = 1:15;
    sdfK            = 1:15;
    extraK          = 0:20;
    nWvMembers      = 5; nWvMembers = nWvMembers*length(monk)*length(areas);
    nFactMembers    = 3; nFactMembers = nFactMembers*length(monk)*length(areas);
    nRespMembers    = 5; nRespMembers = nRespMembers*length(monk)*length(areas);
    compType        = 'all';
    gapRandReps     = 20;
    minGapConsec    = 1;
    gapPerc         = 0;
    gapDim          = 2;
    pcaExplCrit     = 95;
    
    % Waveform clustering parameters
    smoothWaves     = 0;    % Waves smoothed when brought in... see below
    alignWaves      = 0;
    wvRange         = [-250, 500];
    wvClustType     = 'pca';
    waveSimType     = 'corr';
    wvNormType      = 'scale';
    kTypeWv         = 'correlation'; % Use "correlation" or "sqeuclidean"
    pcaTypeWv       = 'agglom';
    pcaAgglomWv     = 'euc';
    pcaKReps        = 50;
    wvGapType       = 'rand';
    wvRandType      = 'pca';
    pcaWaveGapType  = 'rand';
    pcaWaveRandType = 'gauss';
    flipPos         = 0;
    selectK         = 0;
    selWvK          = 0;
    addWave         = 1;
    
    % Response-based clustering params
    respClustType   = 'agglom';
    visSimType      = 'corr';
    movSimType      = 'corr';
    respSimType     = 'euc';
    respNormType    = 'none';
    kTypeVis        = 'correlation';
    kTypeMov        = 'correlation';
    kTypeResp       = 'sqeuclidean';
    doRF            = 0;
    zDim            = 2;
    zType           = 'baseline';
    vWind           = [50:150];
    mWind           = [-100:0];
    visWind         = -200:300; 
    movWind         = -300:200;
    reRandResp      = 1;
    nRepsResp       = 10;
    addResp         = 0;
    selectRespK     = 0;
    gapRespK        = 1;
    visGapType     = 'rand';
    movGapType     = 'rand';
    respGapType     = 'rand';
    visRandType    = 'pca';
    movRandType    = 'pca';
    respRandType    = 'unif';
    
    % Factor-based clustering params
    cutFacts        = {'Pos','Amp','Width','Range','WdSuccess'};
%     cutFacts        = {'Amp','Range','WdSuccess'};
    factSimType     = 'prod';
    factClustType   = 'agglom';
    factNormType    = 'z';
    kTypeFact       = 'sqeuclidean';
    pcaTypeFact     = 'agglom';
    pcaAgglomFact   = 'prod';
    cutCorrFacts    = 1;
    factCorrCutoff  = .1;
    selectF         = 1;
    factGapType     = 'rand';
    factRandType    = 'unif';
    pcaFactGapType  = 'rand';
    pcaFactRandType = 'unif';
    
    % Refinement and self-checking params
    refine          = 1;
    refType         = 'same';
    nReps           = 50;
    nRepsF          = 100;
    
    % Output/Plotting constants
    catNames        = {'WvClust','VisClust','MovClust'};
    printProg       = 1;
    nPrint          = 200;
    plotWaves       = 1;
    plotFacts       = 1;
    plotGroups      = 1;
    plotDend        = 0;
    plotIndivWaves  = 1;
    plotIndivClusts = 1;
    plotWaveIDs     = 0;
    plotWaveFacts   = 0;
    plotRespClusts  = 0;
    testFunFacts    = 0;
    reGroup         = 1;
    plotF           = 0;
    blDim           = 2;
    xlWrite         = 0;
    saveMat         = 1;
    startColor      = [206,153,2]./255;
    midColor        = [120,146,104]./255;
    endColor        = [96,74,123]./255;
    
    breakOutVars = {'MnRate','CV2'}; % Options: MnRate, Fano, CV2
    breakOutTitle = {'Mean Firing Rate (Hz)','CV2'};
    nBreakOutBins = 20;

    % Other constants
%     monk = {'Gauss'};
%     monk = {'Helmholtz'};
    alignEvents = {'StimOnset','SRT'}; goAdj = {'StimOnset','GoCue'};
    xlFile = 'klDataBookKeeping_mg.xlsx';
    sdfWind = {[-200, 1000],[-700, 500]};
    clustAx = [.1 .45 .8 .5];
    sdfAx   = {[.1 .1 .35 .25],[.55 .1 .35 .25]};
    wvTimes         = ((1:32).*25)-(9*25); 
    rawTimes        = wvTimes;
    splTimes        = ((1:.1:32).*25)-(9*25);
    
    areaStr = '';
    for ia = 1:length(areas),
        areaStr = cat(2,areaStr,areas{ia});
    end

    % Load column map
    makeColLookup
    load xlCols
    

    %% Load in the waveforms
    loadTic = tic;
    load masterSorts
    
    fprintf('Extracting waveforms and measured characteristics...');
    if isempty(monkID), reextract = 1; end
    if reextract || ~exist('allWaves','var') || isempty(allWaves),
        allWaveCell = {}; monkID = []; catFacts = []; catCrit = []; funFacts = []; allTypes = {}; allSDF = {}; allTimes = {}; rawSDF = {}; rawNormVals = []; allWvTcell = {}; allAreas = {};
        for im = 1:length(monk),
            [excelNum,~,excelAll] = xlsread(xlFile,monk{im});

            % Get column map for spike measures
            wdCol       = col.spkShape;
            tfrCol      = find(strcmpi(excelAll(4,:),'TfR'),1);
            ampCol      = find(strcmpi(excelAll(4,:),'amplitude'),1);
            relAmpCol   = ampCol + 1;
            mnRateCol   = col.spkStartB + (find(strcmpi(col.spkNamesB,'meanRate'),1)-1);
            fanoCol     = col.spkStartB + (find(strcmpi(col.spkNamesB,'fano'),1)-1);
            cvCol       = col.spkStartB + (find(strcmpi(col.spkNamesB,'cv'),1)-1);
            cv2Col      = col.spkStartB + (find(strcmpi(col.spkNamesB,'cv2'),1)-1);
            lvCol       = col.spkStartB + (find(strcmpi(col.spkNamesB,'lv'),1)-1);
            lvrCol      = col.spkStartB + (find(strcmpi(col.spkNamesB,'lvr'),1)-1);
            mnISICol    = col.spkStartB + (find(strcmpi(col.spkNamesB,'meanISI'),1)-1);
            isoScoreCol = find(strcmpi(excelAll(4,:),'IsoScore'),1);
            fpScoreCol  = isoScoreCol+2;
            snrCol      = col.SNR;
            sigCol      = find(strcmpi(excelAll(4,:),'isSig'),1);
            goodCol     = find(strcmpi(excelAll(4,:),'good'),1);
            
            % Area criterion
            if ismember(areas,'All'),
                 myCrit = logical([zeros(4,1);ones((size(excelAll,1)-4),1)]);
            else
                myCrit = zeros(size(excelAll,1),1);
                for ia = 1:length(areas),
                    myCrit = myCrit | strcmp(excelAll(:,col.area),areas{ia});
                end
            end
            
            if doCuts,
                % Width criterion
                if cutWideNoise,
                    myCrit = myCrit & excelNum(:,wdCol+1) == 1;
                end
                myCrit = myCrit & excelNum(:,wdCol) >= minWidth;

                % ISI criterion
                myCrit = myCrit & (excelNum(:,mnISICol+2) < isiPerc | (excelNum(:,sigCol) == 1));

                % Isolation Score criterion
                myCrit = myCrit & (excelNum(:,isoScoreCol) >= isoMin | (excelNum(:,sigCol) == 1));
                myCrit = myCrit & (excelNum(:,fpScoreCol) <= fpMax | (excelNum(:,sigCol) == 1));
                myCrit = myCrit & excelNum(:,snrCol) >= snrMin;
                
                if sortedOnly,
                    myCrit = myCrit & (excelNum(:,sigCol) == 1) & (excelNum(:,goodCol) == 1);
                end
            end
            
            critInds    = find(myCrit);
            fprintf('\nDiscovered %d valid units for monkey %s...',length(critInds),monk{im}(1:2));
            
            % Get waveforms and SDF
            [monkSDF, monkTimes, monkWaves, monkWvTimes, normMat] = klGetFilev4(critInds,'-m',monk{im},'-t',task,'-rf',doRF,'-z',zDim,'ztype',zType);
            
            % Normalize SDF
            switch respNormType
                case 'none'
                    normSDF = monkSDF;
                case 'max'
                    normSDF = cellfun(@(sdf,blVals,maxVals) (sdf-repmat(blVals(:,1),1,size(sdf,2)))./repmat(maxVals,1,size(sdf,2)), monkSDF, [mat2cell(normMat,ones(size(monkSDF,1),1)),mat2cell(normMat,ones(size(monkSDF,1),1)),mat2cell(normMat,ones(size(monkSDF,1),1))], [mat2cell(normMat(:,3),ones(size(monkSDF,1),1)),mat2cell(normMat(:,4),ones(size(monkSDF,1),1)),mat2cell(normMat(:,3),ones(size(monkSDF,1),1))], 'UniformOutput',0);
                case 'z'
                    normSDF = cellfun(@(sdf,blVals) (sdf-repmat(blVals(:,1),1,size(sdf,2)))./repmat(blVals(:,2),1,size(sdf,2)), monkSDF, [mat2cell(normMat,ones(size(monkSDF,1),1)),mat2cell(normMat,ones(size(monkSDF,1),1)),mat2cell(normMat,ones(size(monkSDF,1),1))],'UniformOutput',0);
                case {'base','bl'}
                    normSDF = cellfun(@(sdf,blMeans) (sdf-repmat(blMeans,1,size(sdf,2))), monkSDF, [mat2cell(normMat(:,1),ones(size(monkSDF,1),1)),mat2cell(normMat(:,1),ones(size(monkSDF,1),1))],'UniformOutput',0);
            end
            rawSDF      = cat(1,rawSDF,monkSDF);
            rawNormVals = cat(1,rawNormVals,normMat);
            allSDF      = cat(1,allSDF,normSDF);
            allTimes    = cat(1,allTimes,monkTimes);
            allWaveCell = cat(1,allWaveCell,monkWaves);
            allWvTcell  = cat(1,allWvTcell,monkWvTimes);

            % Get spike measures from excel file
            rawWd       = excelNum(myCrit,wdCol);
            rawWdSucc   = excelNum(myCrit,wdCol+1);
            rawTFR      = excelNum(myCrit,tfrCol);
            rawAmp      = excelNum(myCrit,ampCol);
            rawRange    = excelNum(myCrit,relAmpCol);
            rawRate     = excelNum(myCrit,mnRateCol);
            rawFano     = excelNum(myCrit,fanoCol);
            rawCV       = excelNum(myCrit,cvCol);
            rawCV2      = excelNum(myCrit,cv2Col);
            rawLV       = excelNum(myCrit,lvCol);
            rawLVR      = excelNum(myCrit,lvrCol);
            rawISI      = excelNum(myCrit,mnISICol);
            rawPos      = rawAmp < 0;
            rawAmp      = abs(rawAmp);

            factMat = [rawPos,rawWd,rawWdSucc,rawAmp,rawRange,rawRate,rawFano,rawCV,rawCV2,rawLV,rawISI];
            factNames = {'Pos','Width','WdSuccess','Amp','Range','MnRate','Fano','CV','CV2','LV','ISI'};
            
            % Get functional properties from excel file
            types       = excelAll(myCrit,col.typeAlt);
            visLat      = excelNum(myCrit,col.visLat);
            movLat      = excelNum(myCrit,col.movLat);
%             visTST      = excelNum(myCrit,col.vTST);
%             movTST      = excelNum(myCrit,col.mTST);
            visTuneDir  = excelNum(myCrit,col.vTune);
            visTuneWd   = excelNum(myCrit,col.vTune+1);
            movTuneDir  = excelNum(myCrit,col.mTune);
            movTuneWd   = excelNum(myCrit,col.mTune+1);


            monkID                  = cat(1,monkID,ones(size(factMat,1),1).*im);
            catFacts                = cat(1,catFacts,factMat);
            catCrit                 = cat(1,catCrit,critInds);
%             funFacts                = cat(1,funFacts,[visLat,movLat,visTST,movTST,visTuneDir,visTuneWd,movTuneDir,movTuneWd]);
%             funNames                = {'VisLat','MovLat','VisTST','MovTST','VisTuneDir','VisTuneWidth','MovTuneDir','MovTuneWidth'};
            allTypes                = cat(1,allTypes,types);
            allAreas                = cat(1,allAreas,excelAll(myCrit,col.area));
        end
        
    else
        allWaves = rawWaves;
        factMat = [rawPos,rawWd,rawWdSucc,rawAmp,rawRange,rawRate,rawFano,rawCV,rawCV2,rawLV,rawISI];
        factNames = {'Pos','Width','WdSuccess','Amp','Range','MnRate','Fano','CV','CV2','LV','ISI'};
        switch respNormType
            case 'none'
                allSDF = rawSDF;
            case 'max'
                allSDF = cellfun(@(sdf,blVals,maxVals) (sdf-repmat(blVals(:,1),1,size(sdf,2)))./repmat(maxVals,1,size(sdf,2)), rawSDF, [mat2cell(rawNormVals,ones(size(rawSDF,1),1)),mat2cell(rawNormVals,ones(size(rawSDF,1),1)),mat2cell(rawNormVals,ones(size(rawSDF,1),1))], [mat2cell(rawNormVals(:,3),ones(size(rawSDF,1),1)),mat2cell(rawNormVals(:,4),ones(size(rawSDF,1),1)),mat2cell(rawNormVals(:,3),ones(size(rawSDF,1),1))], 'UniformOutput',0);
            case 'z'
                allSDF = cellfun(@(sdf,blVals) (sdf-repmat(blVals(:,1),1,size(sdf,2)))./repmat(blVals(:,2),1,size(sdf,2)), rawSDF, [mat2cell(rawNormVals,ones(size(rawSDF,1),1)),mat2cell(rawNormVals,ones(size(rawSDF,1),1)),mat2cell(rawNormVals,ones(size(rawSDF,1),1))],'UniformOutput',0);
            case {'base','bl'}
                allSDF = cellfun(@(sdf,blMeans) (sdf-repmat(blMeans,1,size(sdf,2))), rawSDF, [mat2cell(rawNormVals(:,1),ones(size(rawSDF,1),1)),mat2cell(rawNormVals(:,1),ones(size(rawSDF,1),1))],'UniformOutput',0);
        end
    end
    fprintf('Waveforms extracted in %s\n\n',printTiming(loadTic));
    keepFacts = ~ismember(factNames,cutFacts);
        
    preProcTic = tic;
    fprintf('Beginning pre-processing\n');
    %% Smooth and align waveforms
    if smoothWaves
        fprintf('\tSmoothing waves...');
        if size(allWaves,2) == 32,
            for ii = 1:size(allWaves,1),
                splWaves(ii,:) = spline(1:32,allWaves(ii,:),1:.1:32);
            end
            allWaves = splWaves;
        end
        wvTimes = splTimes;  
        fprintf('Done!\n');
    end

    % Align on derivative
    if alignWaves,
        fprintf('\tAligning Waves...');
        [allWaves, wvTimes] = klTroughAlign(allWaves,wvTimes,0);
        fprintf('Done!\n');
    end
    
    if flipPos,
        allWaves(catFacts(:,1) == 1,:) = -allWaves(catFacts(:,1) == 1,:);
    end
    
% Waves are already smoothed and aligned when brought in... now bring the
% cell into alignment

    wvZero = cellfun(@(times) find(times == 0,1),allWvTcell);
    allWaves = klAlignv2(allWaveCell,wvZero);
    wvTimes = nanmean(klAlignv2(allWvTcell,wvZero),1);
    rawWaves = allWaves;
    
    goodWvInd = wvTimes >= wvRange(1) & wvTimes <= wvRange(2);
    allWaves = allWaves(:,goodWvInd);
    wvTimes = wvTimes(goodWvInd);
    
    %% Normalize waveforms
    switch wvNormType
        case 'zero'
            normWaves = allWaves - repmat(rawWaves(:,1),1,size(allWaves,2));
        case {'norm'}
            minWv = min(allWaves,[],2);
            maxWv = max(allWaves,[],2);
            normWaves = (allWaves-repmat(minWv,1,size(allWaves,2)))./repmat(maxWv-minWv,1,size(allWaves,2));
        case {'scale'},
            maxWv = max(abs(allWaves),[],2);
            normWaves = allWaves./repmat(maxWv,1,size(allWaves,2));
        case {'max'}
            maxWv = max(abs(allWaves-repmat(rawWaves(:,1),1,size(allWaves,2))),[],2);
            normWaves = (allWaves-repmat(rawWaves(:,1),1,size(allWaves,2)))./repmat(maxWv,1,size(allWaves,2));
        case {'mean'}
            normWaves = allWaves - repmat(nanmean(allWaves,2),1,size(allWaves,2));
        otherwise
            normWaves = allWaves;
    end
    
%     if size(wvTimes,2) ~= size(allWaves,2),
%         wvTimes = splTimes;
%     end

    %% Normalize Factor matrix
    switch factNormType
        case {'-z','z'},
            mnFact   = nanmean(catFacts,1);
            stdFact  = nanstd(catFacts,[],1);
            normFacts = (catFacts - repmat(mnFact,size(catFacts,1),1))./repmat(stdFact,size(catFacts,1),1);
        case {'scale','rescale'}
            minFacts = min(catFacts,[],1);
            maxFacts = max(catFacts,[],1);
            normFacts = (catFacts - repmat(minFacts,size(catFacts,1),1))./repmat(maxFacts-minFacts,size(catFacts,1),1);
        otherwise
            normFacts = catFacts;
    end
    
    allFacts = normFacts(:,keepFacts);
    allNames = factNames(keepFacts);
    
    %% Align SDFs
    if reAlignSDF,
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
    end
    fprintf('Pre-Processing Completed in %s\n\n',printTiming(preProcTic));
    
    %% Cluster waveforms
    startTic = tic;
    fprintf('\n*** Clustering Waveforms ***\n');
    if reClust || ~exist('waveClustIDs','var') || isempty(waveClustIDs),
        switch wvClustType
            case 'agglom'
                [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(normWaves,'-k',k,'-r',refine,'-p',printProg,'-t',waveSimType,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint);
                wvSimInd = klGetF(normWaves,waveClustIDs,'-t',waveSimType,'-s',simMat);
            case 'kmeans'
                warning off
                for ir = 1:size(allWaves,1),
                    firstVal(ir) = find(~isnan(normWaves(ir,:)),1);
                    lastVal(ir) = find(~isnan(normWaves(ir,:)),1,'last');
                end
                clpWvs = normWaves(:,max(firstVal):min(lastVal));
                for ik = 1:length(k),
                    for ir = 1:nReps,
                        [icd(:,ir)] = kmeans(clpWvs,k(ik),'Distance',kTypeWv);
                    end
                    waveClustIDs(:,ik) = klModeKmeans(icd);
                    clear icd;
                end
                simMat = [];
                wvSimInd = klGetF(normWaves,waveClustIDs,'-t',waveSimType,'-s',simMat);
            case 'pca'
                warning off
                for ir = 1:size(normWaves,1),
                    firstVal(ir) = find(~isnan(normWaves(ir,:)),1);
                    lastVal(ir) = find(~isnan(normWaves(ir,:)),1,'last');
                end
                clpWvs = normWaves(:,max(firstVal):min(lastVal));
                [~,scores,~,~,expl] = pca(clpWvs);
                cumExpl = cumsum(expl);
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
                        [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(scores,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint,'-w',expl');
                    otherwise
                        [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(scores,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',pcaRefType,'kType',kTypeWv,'nprint',nPrint,'-w',expl');
                end
                wvSimInd = klGetF(scores(:,1:nDims),waveClustIDs,'-t',pcaAgglomWv,'-s',simMat);
%                 [wvK, gapVect, gapErr] = klGetGapv10(normWaves,waveClustIDs,'-t',waveSimType);
            case 'lpp'
                warning off
                for ir = 1:size(allWaves,1),
                    firstVal(ir) = find(~isnan(normWaves(ir,:)),1);
                    lastVal(ir) = find(~isnan(normWaves(ir,:)),1,'last');
                end
                clpWvs = normWaves(:,max(firstVal):min(lastVal));
                [outVals,~,eigVals] = klLPPv1(clpWvs);
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
                        [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(outVals,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint,'-w',eigVals');
                    otherwise
                        [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(outVals,'-k',k,'-r',refine,'-p',printProg,'-t',pcaAgglomWv,'-n',nWvMembers,'refType',pcaRefType,'kType',kTypeWv,'nprint',nPrint,'-w',eigVals');
                end
                wvSimInd = klGetF(scores(:,1:nDims),waveClustIDs,'-t',pcaAgglomWv,'-s',simMat);
%                 [wvK, gapVect, gapErr] = klGetGapv10(normWaves,waveClustIDs,'-t',waveSimType);
            otherwise
                [waveClustIDs, wvClusts, simMat, wvLink] = klAgglomv10(normWaves,'-k',k,'-r',refine,'-p',printProg,'-t',waveSimType,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint);
                wvSimInd = klGetF(normWaves,waveClustIDs,'-t',waveSimType,'-s',simMat);
        end
        % Clear out columns of all-nan
        while sum(isnan(waveClustIDs(:,end))) == size(waveClustIDs,1),
            waveClustIDs = waveClustIDs(:,1:(end-1));
        end
        
    end
    fprintf('\n*** Waveforms Clustered in %s  ***\n',printTiming(startTic));
    
    if strcmp(wvClustType,'pca'),
        [wvK, gapVect, gapErr] = klGetGapv10(scores,waveClustIDs,'-t',pcaAgglomWv,'-n',nFactMembers,'-r',pcaWaveGapType,'randval',pcaWaveRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim,'-w',expl');
    elseif strcmp(wvClustType,'lpp'),
        [wvK, gapVect, gapErr] = klGetGapv10(outVals,waveClustIDs,'-t',pcaAgglomWv,'-n',nFactMembers,'-r',pcaWaveGapType,'randval',pcaWaveRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim,'-w',eigVals');
    else
        [wvK, gapVect, gapErr] = klGetGapv10(normWaves,waveClustIDs,'-t',waveSimType,'-n',nWvMembers,'-r',wvGapType,'randval',wvRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim);
    end

    if length(wvK) > 1,
        fprintf('Found multimodal gap statistics...\n');
        keyboard
    end
    
    if addWave,
        allNames = cat(2,'WvClust',allNames);
        catVars  = ismember(allNames,catNames);
        allFacts = [waveClustIDs(:,k==wvK),allFacts];
    end
    
    %% Cluster SDF
    startTic = tic;
    fprintf('\n*** Clustering SDFs ***\n');f
    if reClustResp %|| ~exist('visClustIDs','var') || isempty(visClustIDs),
        [respMat, respNames, respSlope] = klParseEpochsv2(visAlign, visTimes, movAlign, movTimes,'-s',1);
        switch respClustType
            case 'agglom'
                [visClustIDs, visClusts, visSimMat, visLink] = klAgglomv10(visAlign(:,ismember(nanmean(visTimes,1),vWind)),'-k',sdfK,'-r',refine,'-p',printProg,'-t',visSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeVis,'nprint',nPrint);
                [movClustIDs, movClusts, movSimMat, movLink] = klAgglomv10(movAlign(:,ismember(nanmean(movTimes,1),mWind)),'-k',sdfK,'-r',refine,'-p',printProg,'-t',movSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeMov,'nprint',nPrint);
                [respClustIDs, respClusts, respSimMat, respLink] = klAgglomv10(respMat,'-k',sdfK,'-r',refine,'-p',printProg,'-t',respSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeResp,'nprint',nPrint);
                while sum(isnan(visClustIDs(:,end))) == size(visClustIDs,1),
                    visClustIDs = visClustIDs(:,1:(end-1));
                end
                while sum(isnan(movClustIDs(:,end))) == size(movClustIDs,1),
                    movClustIDs = movClustIDs(:,1:(end-1));
                end
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
            otherwise
                [visClustIDs, visClusts, visSimMat, visLink] = klAgglomv10(visAlign(:,ismember(nanmean(visTimes,1),vWind)),'-k',sdfK,'-r',refine,'-p',printProg,'-t',respSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeResp,'nprint',nPrint);
                [movClustIDs, movClusts, movSimMat, movLink] = klAgglomv10(movAlign(:,ismember(nanmean(movTimes,1),mWind)),'-k',sdfK,'-r',refine,'-p',printProg,'-t',respSimType,'-n',nRespMembers,'refType',refType,'kType',kTypeResp,'nprint',nPrint);
        end
        [visSimF,visSimR] = klGetF(visAlign(:,ismember(nanmean(visTimes,1),visWind)),visClustIDs,'-t',visSimType,'-s',visSimMat);
        [movSimF,movSimR] = klGetF(movAlign(:,ismember(nanmean(movTimes,1),movWind)),movClustIDs,'-t',movSimType,'-s',movSimMat);
        [respSimF, respSimR] = klGetF(respMat,respClustIDs,'-t',respSimType,'-s',respSimMat);
        
        [visK,visGap,visGapErr] = klGetGapv10(visAlign(:,ismember(nanmean(visTimes,1),visWind)),visClustIDs,'-t',visSimType,'-r',visGapType,'randval',visRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim);
        [movK,movGap,movGapErr] = klGetGapv10(movAlign(:,ismember(nanmean(movTimes,1),movWind)),movClustIDs,'-t',movSimType,'-r',movGapType,'randval',movRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim);
        [respK,respGap,respGapErr] = klGetGapv10(respMat,respClustIDs,'-t',respSimType,'-r',respGapType,'randval',respRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim);
        
    end
    
    fprintf('\n*** SDF Clustered in %s  ***\n',printTiming(startTic));
    
    if addResp,
        allNames = cat(2,'VisClust',allNames);
        allNames = cat(2,'MovClust',allNames);
        allFacts = [visClustIDs(:,visK),allFacts];
        allFacts = [movClustIDs(:,movK),allFacts];
    end
    
    if cutCorrFacts,
        allFactsPre = allFacts;
        allNamesPre = allNames;
        sigCorrCut = [];
        keepCorr = [];
        startCol = find(strcmp(allNames,'MnRate'));
        for iFact1 = 1:(length(allNames)+1-startCol),
            for iFact2 = (iFact1+1):(length(allNames)+1-startCol),
                clear fact1 fact2
                fact1 = allFacts(:,startCol+(iFact1-1));
                fact2 = allFacts(:,startCol+iFact2-1);
%                 [r,p] = corr(fact1(~isnan(fact1) & ~isnan(fact2)),fact1(~isnan(fact1) & ~isnan(fact2)));
                [r,p] = corr(fact1,fact2);
                if p < factCorrCutoff && ~ismember(startCol-1+iFact2,sigCorrCut),
                    sigCorrCut(length(sigCorrCut)+1) = startCol+iFact2-1;
                    keepCorr(length(keepCorr)+1) = startCol+(iFact1-1);
                end
            end
        end
        allNames = allNames(~ismember(1:size(allFacts,2),sigCorrCut));
        allFacts = allFacts(:,~ismember(1:size(allFacts,2),sigCorrCut));
    end
    catVars  = ismember(allNames,catNames);
    
    if reClust,
        fprintf('Factor-based clustering with %d wave clusters...\n',wvK);
            
        % Do the agglomeration with new factors included
        switch factClustType
            case 'agglom'
                [factClustsOut, factClustMembs, factSim, factLink] = klAgglomv10(allFacts,'-k',(wvK+extraK)','-r',refine,'-p',printProg,'-t',factSimType,'-c',compType,'-cv',catVars,'nprint',nPrint,'-n',nFactMembers);
            case 'kmeans'
                for iik = 1:length(extraK),
                    for ir = 1:nReps,
                        [icd(:,ir)] = kmeans(allFacts,wvK+extraK(iik),'Distance',kTypeFact);
                    end
                    factClustsOut(:,iik) = klModeKmeans(icd);
                    clear icd;
                end
            case 'pca',
                goodRows = ones(size(allFacts,1),1);
                for i = 1:size(allFacts,2),
                    goodRows = goodRows & ~isnan(allFacts(:,i));
                end
                [coeffs,scores,~,~,expl] = pca(allFacts(goodRows,~catVars));
                cumExpl = cumsum(expl);
                nDims = find(cumExpl > pcaExplCrit,1);
                allScores = nan(size(allFacts,1),size(scores,2));
                allScores(goodRows,:) = scores;
                pcaFacts(:,1:sum(catVars)) = allFacts(:,catVars);
                pcaFacts(:,(sum(catVars)+1):(sum(catVars)+nDims)) = allScores(:,1:nDims);
                
                switch pcaTypeFact,
                    case 'agglom'
                        [factClustsOut,factClustMembs,factSim,factLink] = klAgglomv10(pcaFacts,'-k',(wvK+extraK)','-r',refine,'-p',printProg,'-t',pcaAgglomFact,'-c',compType,'-cv',catVars,'nprint',nPrint,'-n',nFactMembers);
                end
                
            otherwise
                [factClustsOut, factClustMembs, factSim, factLink] = klAgglomv10(allFacts,'-k',(wvK+extraK)','-r',refine,'-p',printProg,'-t',factSimType,'-c',compType,'-cv',catVars,'nprint',nPrint,'-n',nFactMembers);
        end
        [clustDist, clustR, clustBet, clustWit] = klGetF(allFacts,factClustsOut,'-t',factSimType,'cat',catVars,'-s',factSim);
        
        % Clear out columns of all-nan
        while sum(isnan(factClustsOut(:,end))) == size(factClustsOut,1),
            factClustsOut = factClustsOut(:,1:(end-1));
        end
        
        % Randomly assign cluster IDs and check F again
        clear randDist randR randBet randWit
        for ii = 1:nRepsF
            if mod(ii,50) == 0, fprintf('%d loops done\n',ii); end
            randInt = randperm(size(factClustsOut,1));
            [randDist(ii,:), randR(ii,:), randBet(ii,:), randWit(ii,:)] = klGetF(allFacts,factClustsOut(randInt,:),'-t',factSimType,'cat',catVars,'-s',factSim(randInt,randInt));  
        end
        
        if strcmp(factClustType,'pca'),
            [factK, factGap, factErr] = klGetGapv10(pcaFacts,factClustsOut,'-t',pcaAgglomFact,'-n',nFactMembers,'cat',catVars,'-r',pcaFactGapType,'randval',pcaFactRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim);
        else
            [factK, factGap, factErr] = klGetGapv10(allFacts,factClustsOut,'-t',factSimType,'-n',nFactMembers,'cat',catVars,'-r',factGapType,'randval',factRandType,'randreps',gapRandReps,'-c',minGapConsec,'perc',gapPerc,'dim',gapDim);
        end
        factClustIDs = factClustsOut(:,(wvK+extraK) == factK);
    end
end

reArrange = nan(size(factClustIDs));
nDone = 1;
for iw = 1:wvK,
    thisGroup = unique(factClustIDs(waveClustIDs(:,k==wvK) == iw));
    for ig = 1:length(thisGroup),
        reArrange(factClustIDs == thisGroup(ig)) = nDone;
        nDone = nDone+1;

    end
end
factClustIDs = reArrange;
    
% clustColors = jet(factK);
% respColors = jet(max([visK,movK]));
wvClustColors = klCMapv2(startColor,midColor,endColor,wvK);
clustColors = klCMapv2(startColor,midColor,endColor,factK);
respColors = klCMapv2(startColor,midColor,endColor,max([visK,movK,respK]));
if saveMat
    save('./clustVars_new.mat');
end

if flipPos,
    normWaves(catFacts(:,1) == 1,:) = -normWaves(catFacts(:,1) == 1,:);
end
    
%% Plot stuff
% Plot waveforms and mean waveforms for the clusters
if plotWaves
    figure(200);
    myWaveIDs = waveClustIDs(:,k == wvK);
    uWave = unique(myWaveIDs(~isnan(myWaveIDs)));
    mnWave = nan(length(uWave),size(normWaves,2));
    for iw = 1:length(uWave),
        mnWave(iw,:) = nanmean(normWaves(myWaveIDs == uWave(iw),:),1);
        wvLegStr{iw} = sprintf('Clus. %d (n=%d)',iw,sum(myWaveIDs==uWave(iw)));
        plot(wvTimes,mnWave(iw,:),'color',wvClustColors(iw,:),'linewidth',3); hold on;
    end
    set(gca,'fontsize',16,'tickdir','out','XLim',[-200, 500]);
    legend(wvLegStr);
    xlabel('Time (us)','fontsize',20);
    ylabel('Relative Voltage','fontsize',20);
    title('Mean Within Cluster Waveforms','fontsize',24);

end
        
if plotFacts,
    for ik = 1:size(factClustIDs,2),
        clear clustMeans
        clustInds = unique(factClustIDs(:,ik));
        clustInds = clustInds(~isnan(clustInds));
        endLegStr = {};
        for ic = 1:length(clustInds),
            waveInds = factClustIDs(:,ik) == clustInds(ic);
            clustMeans(ic,:) = nanmean(normWaves(waveInds,:),1);
            endLegStr{end+1} = sprintf('Clus. %d (n=%d)',ic,sum(waveInds));
            if plotIndivClusts
                figure(ic);
                if plotGroups, axes('position',clustAx); end
                if plotIndivWaves,
                    plot(wvTimes,normWaves(waveInds,:)');
                    hold on;
                    plot(wvTimes,clustMeans(ic,:),'k','linewidth',3);
                else
                    pltMeanStd(wvTimes,clustMeans(ic,:),nanstd(normWaves(waveInds,:),[],1)./sqrt(sum(isfinite(normWaves(waveInds,:)),1)),'k');
                end
                t=title(sprintf('Wave Cluster %d (n=%d)',ic,sum(waveInds)));
            end
        end

        figure(100);
        for ic = 1:size(clustMeans,1),
            plot(wvTimes,clustMeans(ic,:),'color',clustColors(ic,:),'linewidth',3);
            hold on;
        end
        set(gca,'XTick',0:200:800,'fontsize',16,'tickdir','out','XLim',[-200, 500]);
        xlabel('Time (us)','fontsize',20);
        ylabel('Relative Voltage','fontsize',20);
        title('Mean Within Cluster Waveforms','fontsize',24);
        legend(endLegStr,'fontsize',16);

        if plotDend
            figure(1000+ik);
            dendrogram(wvLink,0);
        end
    end
end

if plotGroups,
    for ik = 1:1%size(factClustIDs,2),
        clustInds = unique(factClustIDs(:,ik));
        legStr = {};
        for ic = 1:length(clustInds),
            fprintf('Plotting Cluster %d (%d members)\n',ic,sum(factClustIDs(:,ik) == clustInds(ic)));
            
            % Plot visual responses
            figure(ic); 
            axes('position',sdfAx{1});
            myVis = visAlign(factClustIDs(:,ik)==clustInds(ic),ismember(nanmean(visTimes,1),visWind));
            pltMeanStd(nanmean(visTimes(:,ismember(nanmean(visTimes,1),visWind)),1),nanmean(myVis,1),nanstd(myVis,[],1)./sqrt(sum(isfinite(myVis),1)),'k');
            vline(0);
            set(gca,'XLim',[visWind(1),visWind(end)],'tickdir','out');
            xlabel('Time from Stimulus On'); ylabel('SDF');
            
            axes('position',sdfAx{2});
            myMov = movAlign(factClustIDs(:,ik)==clustInds(ic),ismember(movTimes,movWind));
            pltMeanStd(movTimes(:,ismember(movTimes,movWind)),nanmean(myMov,1),nanstd(myMov,[],1)./sqrt(sum(isfinite(myMov),1)),'k');
            vline(0);
            set(gca,'XLim',[movWind(1),movWind(end)],'tickdir','out');
            xlabel('Time from Saccade'); ylabel('SDF');
            
            figure(300+ik);
            subplot(1,2,1);
            plot(nanmean(visTimes(:,ismember(nanmean(visTimes,1),visWind)),1),nanmean(myVis,1),'color',clustColors(ic,:)); hold on;
            legStr = cat(2,legStr,sprintf('Cluster %d',ic));
            if ic == length(clustInds),
%                 set(gca,'color','k');
                vvl = vline(0); set(vvl,'color','k');
                xlabel('Time from Stim Onset','fontsize',18);
                ylabel('Normalized SDF','fontsize',18);
                legend(legStr,'color','w');
            end
            subplot(1,2,2);
            plot(nanmean(movTimes(:,ismember(nanmean(movTimes,1),movWind)),1),nanmean(myMov,1),'color',clustColors(ic,:)); hold on;
            if ic == length(clustInds),
%                 set(gca,'color','k');
                mvl = vline(0); set(mvl,'color','k');
                xlabel('Time from SRT','fontsize',18);
                ylabel('Normalized SDF','fontsize',18);
                legend(legStr);
                st = suptitle(sprintf('Mean SDF by Cluster')); set(st,'fontsize',26);
            end
            
        end
    end
end

if plotWaveFacts,
    minVar1 = min(catFacts(:,strcmp(factNames,breakOutVars{1})));
    maxVar1 = max(catFacts(:,strcmp(factNames,breakOutVars{1})));
    minVar2 = min(catFacts(:,strcmp(factNames,breakOutVars{2})));
    maxVar2 = max(catFacts(:,strcmp(factNames,breakOutVars{2})));

    varBins1 = minVar1:(range(catFacts(:,strcmp(factNames,breakOutVars{1})))/nBreakOutBins):maxVar1;
    varBins2   = minVar2:(range(catFacts(:,strcmp(factNames,breakOutVars{2})))/nBreakOutBins):maxVar2;

    for ic = 1:wvK,
         % Open figure and axes
         figure(1000+ic);
         wvAx = axes('position',clustAx);
         varAx1 = axes('position',sdfAx{1});
         varAx2 = axes('position',sdfAx{2});

         % Get appropriate rows
         myRows = waveClustIDs(:,k==wvK)==ic;

         % Get fact clusts
         childFacts = unique(factClustIDs(myRows)); childFacts(isnan(childFacts)) = [];

         clear varHist1 varHist2
         for iic = 1:length(childFacts),
             axes(wvAx);
             plot(wvTimes,normWaves(factClustIDs == childFacts(iic),:),'color',clustColors(childFacts(iic),:));
             hold on;

             axes(varAx1);
             varHist1(iic,:) = hist(catFacts(factClustIDs == childFacts(iic),strcmp(factNames,breakOutVars{1})),varBins1);
             hold on;

             axes(varAx2);
             varHist2(iic,:) = hist(catFacts(factClustIDs == childFacts(iic),strcmp(factNames,breakOutVars{2})),varBins2);
             hold on;
         end

         axes(wvAx);
         xlabel('Time (us)','fontsize',12); ylabel('Normalized Voltage (A.U.)','fontsize',12);
         t=title(sprintf('Child Clusters of Waveform Cluster %d',ic),'fontsize',18);

         axes(varAx1);
         rateBar = bar(varBins1,varHist1',1.2);
         for iic = 1:length(childFacts),
             set(rateBar(iic),'FaceColor',clustColors(childFacts(iic),:));
         end
         xlabel(breakOutTitle{1},'fontsize',10);
         ylabel('Count');

         axes(varAx2);
         cvBar = bar(varBins2,varHist2',1.2);
         for iic = 1:length(childFacts),
             set(cvBar(iic),'FaceColor',clustColors(childFacts(iic),:));
         end
         xlabel(breakOutTitle{2},'fontsize',10);
         ylabel('Count');

    end         
end

if plotRespClusts,
    for ic = 1:visK,
        figure(100+ic);
        axes('position',clustAx);
        plot(wvTimes,normWaves(visClustIDs(:,sdfK == visK) == ic,:));
        hold on;
        plot(wvTimes,nanmean(normWaves(visClustIDs(:,sdfK == visK) == ic,:),1),'color','k','linewidth',2);
        xlabel('Time (us)','fontsize',18); ylabel('Relative Voltage','fontsize',18);
        t=title(sprintf('Vis Cluster %d (n=%d)',ic,sum(visClustIDs(:,sdfK == visK) == ic)));
                
        axes('position',sdfAx{1});
%         plot(visWind,nanmean(visAlign(visClustIDs(:,respK == visK) == ic,ismember(nanmean(visTimes,1),visWind)),1),'color',respColors(ic,:));
        pltMeanStd(visWind,nanmean(visAlign(visClustIDs(:,sdfK == visK) == ic,ismember(nanmean(visTimes,1),visWind)),1),nanstd(visAlign(visClustIDs(:,sdfK == visK) == ic,ismember(nanmean(visTimes,1),visWind)),[],1)./sqrt(sum(isfinite(visAlign(visClustIDs(:,sdfK == visK) == ic,ismember(nanmean(visTimes,1),visWind))),1)),'k');
        vline(0);
        xlabel('Time from StimOn'); ylabel('Normalized SDF');
        
        axes('position',sdfAx{2});
%         plot(movWind,nanmean(movAlign(visClustIDs(:,respK == visK) == ic,ismember(nanmean(movTimes,1),movWind)),1),'color',respColors(ic,:));
        pltMeanStd(movWind,nanmean(movAlign(visClustIDs(:,sdfK == visK) == ic,ismember(nanmean(movTimes,1),movWind)),1),nanstd(movAlign(visClustIDs(:,sdfK == visK) == ic,ismember(nanmean(movTimes,1),movWind)),[],1)./sqrt(sum(isfinite(movAlign(visClustIDs(:,sdfK == visK) == ic,ismember(nanmean(movTimes,1),movWind))),1)),'k');
        vline(0);
        xlabel('Time from SRT'); ylabel('Normalized SDF');
    end
    for ic = 1:movK,
        figure(1000+ic);
        axes('position',clustAx);
        plot(wvTimes,normWaves(movClustIDs(:,sdfK == visK) == ic,:));
        hold on;
        plot(wvTimes,nanmean(normWaves(movClustIDs(:,sdfK == movK) == ic,:),1),'color','k','linewidth',2);
        xlabel('Time (us)','fontsize',18); ylabel('Relative Voltage','fontsize',18);
        t=title(sprintf('Mov Cluster %d (n=%d)',ic,sum(movClustIDs(:,sdfK == movK) == ic)));
                
        axes('position',sdfAx{1});
%         plot(visWind,nanmean(visAlign(visClustIDs(:,respK == visK) == ic,ismember(nanmean(visTimes,1),visWind)),1),'color',respColors(ic,:));
        pltMeanStd(visWind,nanmean(visAlign(movClustIDs(:,sdfK == movK) == ic,ismember(nanmean(visTimes,1),visWind)),1),nanstd(visAlign(movClustIDs(:,sdfK == movK) == ic,ismember(nanmean(visTimes,1),visWind)),[],1)./sqrt(sum(isfinite(visAlign(movClustIDs(:,sdfK == movK) == ic,ismember(nanmean(visTimes,1),visWind))),1)),'k');
        vline(0);
        xlabel('Time from StimOn'); ylabel('Normalized SDF');
        
        axes('position',sdfAx{2});
%         plot(movWind,nanmean(movAlign(visClustIDs(:,respK == visK) == ic,ismember(nanmean(movTimes,1),movWind)),1),'color',respColors(ic,:));
        pltMeanStd(movWind,nanmean(movAlign(movClustIDs(:,sdfK == movK) == ic,ismember(nanmean(movTimes,1),movWind)),1),nanstd(movAlign(movClustIDs(:,sdfK == movK) == ic,ismember(nanmean(movTimes,1),movWind)),[],1)./sqrt(sum(isfinite(movAlign(movClustIDs(:,sdfK == movK) == ic,ismember(nanmean(movTimes,1),movWind))),1)),'k');
        vline(0);
        xlabel('Time from SRT'); ylabel('Normalized SDF');
    end
end
        
if testFunFacts,
    for i = 1:size(funFacts,2),
        [pw(i),~,statsw] = kruskalwallis(funFacts(:,i),waveClustIDs(:,k == wvK),'off');
        [pf(i),~,statsf] = kruskalwallis(funFacts(:,i),factClustIDs,'off');
        if pw(i) < .05,
            figure(2000+i); multcompare(statsw); title(sprintf('%s by WF Cluster: p=%.4f',funNames{i},pw(i)));
        end
        if pf(i) < .05,
            figure(3000+i); multcompare(statsf); title(sprintf('%s by Factor-Based Cluster: p=%.4f',funNames{i},pf(i)));
        end
    end
end
        
if plotWaveIDs
    % Reality check for categorical clustering appropriateness
    for ii = 1:factK
        figure(10000+ii);
        legStr = {};
        for ic = 1:wvK
            if sum(allFacts(:,1) == ic & factClustIDs == ii) > 0,
                plot(wvTimes,normWaves(allFacts(:,1) == ic & factClustIDs == ii,:),'color',clustColors(ic,:));
                legStr = cat(2,legStr,sprintf('WvClust %d (n=%d)',ic,sum(allFacts(:,1) == ic & factClustIDs == ii)));
                hold on;
            end
        end
        legend(legStr);
    end
end

if xlWrite
    lastCol  = find(strcmp(excelAll(4,:),'lastClust'));
    wvHead  = ['wvClust',areaStr];
    factHead = ['factClust',areaStr];
    visHead = ['visClust',areaStr];
    movHead = ['movClust',areaStr];
    wvClustCol = find(strcmp(excelAll(4,:),wvHead));
    factClustCol = find(strcmp(excelAll(4,:),factHead));
    visClustCol = find(strcmp(excelAll(4,:),visHead));
    movClustCol = find(strcmp(excelAll(4,:),movHead));
    
    for im = 1:length(monk),
        % Write the waves
        wvWrite     = nan(max(catCrit(monkID == im)),1);
        wvWrite(catCrit(monkID == im)) = waveClustIDs(monkID == im,k == wvK);
        if ~isempty(wvClustCol),
            xlswrite(xlFile,wvWrite,monk{im},sprintf('%s1',num2abc(wvClustCol)));
            xlswrite(xlFile,{wvHead},monk{im},sprintf('%s4',num2abc(wvClustCol)));
        else
            warning('Unable to find appropriate column. Wave clusters not written');
        end
        
        % Write the visual response
        visWrite     = nan(max(catCrit(monkID == im)),1);
        visWrite(catCrit(monkID == im)) = visClustIDs(monkID == im, sdfK == visK);
        if ~isempty(visClustCol),
            xlswrite(xlFile,visWrite,monk{im},sprintf('%s1',num2abc(visClustCol)));
            xlswrite(xlFile,{visHead},monk{im},sprintf('%s4',num2abc(visClustCol)));
        else
            warning('Unable to find appropriate column. Visual clusters not written');
        end
        
        % Write the movement response
        movWrite     = nan(max(catCrit(monkID == im)),1);
        movWrite(catCrit(monkID == im)) = movClustIDs(monkID == im, sdfK == movK);
        if ~isempty(movClustCol),
            xlswrite(xlFile,movWrite,monk{im},sprintf('%s1',num2abc(movClustCol)));
            xlswrite(xlFile,{movHead},monk{im},sprintf('%s4',num2abc(movClustCol)));
        else
            warning('Unable to find appropriate column. Movement clusters not written');
        end
        
        % Write the factor-based clusters
        factWrite     = nan(max(catCrit(monkID == im)),1);
        factWrite(catCrit(monkID == im)) = factClustIDs(monkID == im);
        xlswrite(xlFile,factWrite,monk{im},sprintf('%s1',num2abc(lastCol)));
        xlswrite(xlFile,{'lastClust'},monk{im},sprintf('%s4',num2abc(lastCol)));
        if ~isempty(factClustCol),
            xlswrite(xlFile,factWrite,monk{im},sprintf('%s1',num2abc(factClustCol)));
            xlswrite(xlFile,{factHead},monk{im},sprintf('%s4',num2abc(factClustCol)));
        else
            warning('Unable to find appropriate column. Factor clusters written to "last" only');
        end
    end
end

fprintf('\nScript completed in %s\n',printTiming(grandTic));