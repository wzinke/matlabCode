function [respK, respClustIDs, respGap, respGapErr, linkMat, linkInd, thisDistOut, respInOut] = klRespClustForMetaclustv2(allSDF,allTimes,varargin)

% global sdfK refine printProg respSimType nRespMembers refType kTypeResp nPrint respClustType pcaExplCrit pcaTypeWv pcaAgglomWv pcaWaveGapType pcaWaveRandType gapRandReps minGapConsec gapPerc gapDim
% load('origClustParams.mat');

%% Get aligned SDFs and Times
if ~any(size(allSDF)==1),
    fprintf('\tAligning SDF and times...');
    visZero = cellfun(@(times) find(times == 0,1),allTimes(:,1));
    movZero = cellfun(@(times) find(times == 0,1),allTimes(:,2));
    visAlign = klAlignv2(allSDF(:,1),visZero);
    visTimes = nanmean(klAlignv2(allTimes(:,1),visZero),1);
    movAlign = klAlignv2(allSDF(:,2),movZero);
    movTimes = nanmean(klAlignv2(allTimes(:,2),movZero),1);
    visAlign(visAlign == inf) = nan;
    movAlign(movAlign == inf) = nan;
    fprintf('Done!\n');
else
    visAlign = allSDF{1}; visTimes = allTimes{1};
    movAlign = allSDF{2}; movTimes = allTimes{2};
end

%% Set Defaults
respNormType = 'none';
respSimType = 'corr';
respInclude = 'combo';
zResp = 1;
sdfK = 1:20;
blWind = -300:-100;
vWind = -200:300;
mWind = -300:200;
nRespMembers = 10;
nSkip = 1;
randReps = 10; 
direction = 'first';
dirTypes = {'first','last'};
doGap = 1;

% Default analysis epochs
preVis = -100:0;
visTrans = 50:100;
visSust = 100:150;
preMov = -50:0;
postMov = 0:50;
nextVis = 50:100;
allEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis};



%% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'norm'},
            respNormType = varargin{varStrInd(iv)+1};
        case {'sim'},
            respSimType = varargin{varStrInd(iv)+1};
        case {'resp'},
            respInclude = varargin{varStrInd(iv)+1};
        case {'z'},
            zResp = varargin{varStrInd(iv)+1};
        case {'blwind'},
            blWind = varargin{varStrInd(iv)+1};
        case {'-n'},
            nRespMembers = varargin{varStrInd(iv)+1};
        case {'-d'},
            direction = varargin{varStrInd(iv)+1};
            if isnumeric(direction),
                direction = dirTypes{direction};
            end
        case {'-g','gap'},
            doGap = varargin{varStrInd(iv)+1};
        case {'-r','reps'},
            randReps = varargin{varStrInd(iv)+1};
    end
end

%% Start the analysis

% Normalize/Scale the SDFs
[normResp{1}, normMov] = klNormResp(visAlign,visTimes,movAlign,movTimes,respNormType,'bl',blWind);

% Cut down the responses to relevant times and to get rid of nans
catNorms = [normResp{1}(:,ismember(visTimes,vWind)),normMov(:,ismember(movTimes,mWind))];
goodRows = sum(isfinite(catNorms),2) == size(catNorms,2);
normResp{1} = normResp{1}(goodRows,:);
normMov = normMov(goodRows,:);

epocSpks = {normVis,normVis,normVis,normMov,normMov,normMov};
epocTms = {visTimes,visTimes,visTimes,movTimes,movTimes,movTimes};

% Smooth the SDFs a little bit in case we want to cluster the whole SDF
kern = klGetKern('type','gauss','width',10);
convVis = conv2(normResp{1}(:,ismember(visTimes,vWind)),kern,'same');
convMov = conv2(normMov(:,ismember(movTimes,mWind)),kern,'same');
convResp = [convVis(:,1:nSkip:end),convMov(:,1:nSkip:end)];
startInds = [1,length(1:nSkip:size(convVis,2))+1];

%% Get the proper input ready
startTic = tic;
fprintf('\n*** Clustering SDFs ***\n');
[respMat, respSlope] = klParseEpochsv3(epocSpks,epocTms,allEpocs);

switch respInclude
    case {'mean','means'}
        respIn = respMat;
    case {'slope','slopes'}
        respIn = respSlope;
    case {'comb','combo'},
        respIn = [respMat,respSlope];
    case {'whole'},
        respIn = convResp;
end
if zResp,
    respIn = (respIn-repmat(nanmean(respIn,1),size(respIn,1),1))./repmat(nanstd(respIn,[],1),size(respIn,1),1);
end

%% Make a distance matrix
thisDist = makeDistMat(respIn,respSimType,zeros(1,size(respIn,2)));
thisDistOut = nan(size(visAlign,1));
thisDistOut(goodRows,goodRows) = thisDist;

%% Cluster the distance matrix
[respClustIDs,~,~,clusMembers,linkMat, linkInd] = klDistMatAgglomv2(thisDist,sdfK,'-n',nRespMembers);
nBad = find(~goodRows);
for i = 1:length(nBad)
    check = linkMat(:,1:2);
    linkMat(check >= nBad(i)) = linkMat(check >= nBad(i))+1;
    linkMat = cat(1,linkMat,[nBad(i),max(check(:))+1,linkMat(end,3)]);
end
nMembs = cellfun(@length,clusMembers);
nGood = sum(nMembs >= nRespMembers,1);
while sum(isnan(respClustIDs(:,end))) == size(respClustIDs,1)
    respClustIDs = respClustIDs(:,1:(end-1));
end

% Calculate W for each value of k (see Tibshirani et al 2001 for gap
% statistic calculations). 
W = nan(1,size(respClustIDs,2));
for ik = 1:size(respClustIDs,2),
    W(ik) = getW(respClustIDs(:,ik),thisDist);
end

if doGap,
    %% Get the gap statistic (See Tibshirani et al 2001)
    wStar = nan(randReps,size(respClustIDs,2));
    for iRand = 1:randReps,
        printStr = sprintf('Working on random repetition %d of %d...',iRand,randReps);
        fprintf(printStr);

        % Get a random set of data according to the mean/SD of the inputs
        randMat = klMakeRand(size(respIn,1),nanmean(respIn,1),nanstd(respIn,1));
        randDist = makeDistMat(randMat,respSimType,zeros(1,size(randMat,2)));
        randIDs = klDistMatAgglomv2(randDist,sdfK,'-n',nRespMembers,'-p',0);
        while sum(isnan(randIDs(:,end))) == size(randIDs,1),
            randIDs = randIDs(:,1:(end-1));
        end
        % Calculate randomized W
        for ik = 1:min([size(wStar,2),size(randIDs,2)]),
            wStar(iRand,ik) = getW(randIDs(:,ik),randDist);
        end
        for ib = 1:length(printStr),
            fprintf('\b');
        end
    end
    fprintf('Calculating on gap...\n');

    % First, calculate the gap itself
    respGap = nanmean(log(wStar)-repmat(log(W),size(wStar,1),1),1);

    % Now we need lBar
    lBar = nanmean(log(wStar),1);

    % Now s(respGapErr)/sd...
    sd = sqrt(nanmean((log(wStar)-repmat(lBar,size(wStar,1),1)).^2,1));
    respGapErr = sd.*sqrt(1+(1/randReps));

    % Now get the "official" k
    gapComp = [respGap(2:end)-respGapErr(2:end),nan];
    respK = find(respGap > gapComp,1,direction);
    if isempty(respK),
        respK = nan;
    end
else
    respGap = nan;
    respGapErr = nan;
    respK = nan;
end

% Set the variables aside...
respClustTmp = respClustIDs;
normVisTmp = normVis; normMovTmp = normMov;
clear respClustIDs normVis normMov;

% Initialize the output variables
respClustIDs = nan(size(catNorms,1),size(respClustTmp,2));
normVis = nan(size(catNorms,1),size(normVisTmp,2));
normMov = nan(size(catNorms,1),size(normMovTmp,2));
respMatOut = nan(size(catNorms,1),size(respMat,2));
respInOut = nan(size(catNorms,1),size(respIn,2));

% Place values in the output variables
respClustIDs(goodRows,:) = respClustTmp;
normVis(goodRows,:) = normVisTmp;
normMov(goodRows,:) = normMovTmp;
respMatOut(goodRows,:) = respMat;
respInOut(goodRows,:) = respIn;

%% And we're done!
fprintf('\n*** SDF Clustered in %s  ***\n',printTiming(startTic));


function distMat = makeDistMat(obs,type,catVars,print)
    if nargin < 4,
        print = 0;
    end
    % Rescale again for product distance
    if ismember(type,{'prod'}),
        maxObs = nanmax(obs,[],1);
        minObs = nanmin(obs,[],1);
        normObs = (obs-repmat(minObs,size(obs,1),1))./repmat((maxObs-minObs),size(obs,1),1);
        normRaw = normObs;
    end

    % Start by computing pair-wise distances
    if print,
        fprintf('Initializing pairwise distances...')
    end
    if ismember(type,{'corr'}),
        % Clean up obs, if "corr"
        goodCols = zeros(1,size(obs,2));
        for ic = 1:size(obs,2),
            goodCols(ic) = sum(~isnan(obs(:,ic))) == size(obs,1);
        end
        obs = obs(:,logical(goodCols));
        distMat = 1-corr(obs',obs');
        for ir = 1:size(obs,1),
            for ic = ir:size(obs,1),
                distMat(ic,ir) = nan;
            end
        end
    elseif ismember(type,{'euc'}),
        distMat = EuDist2(obs(:,~catVars));
        for ir = 1:size(obs,1),
            for ic = ir:size(obs,1),
                distMat(ic,ir) = nan;
            end
        end
    else
        distMat = nan(size(obs,1),size(obs,1));
        for im = 1:size(obs,1),
            for in = (im+1):size(obs,1),
                if ismember(type,{'prod','exprod'}) || doNorm,
                    distMat(im,in) = klGetDistv2(normObs(im,:),normObs(in,:),'-t',type,'cat',catVars);
                else
                    distMat(im,in) = klGetDistv2(obs(im,:),obs(in,:),'-t',type,'cat',catVars);
                end
            end
        end
    end
    if print,
        fprintf('Done!\n');
    end
end

function timeStr = printTiming(ticVal)
    thisRawTime = toc(ticVal);
    thisTimeMin = floor(thisRawTime/60);
    if thisTimeMin > 60, thisTimeHour = floor(thisTimeMin/60); thisTimeMin = mod(thisTimeMin,60); else thisTimeHour = 0; end
    thisTimeSec = mod(thisRawTime,60);

    thisHrStr = num2str(thisTimeHour); while length(thisHrStr) < 2, thisHrStr = ['0',thisHrStr]; end
    thisMinStr = num2str(thisTimeMin); while length(thisMinStr) < 2, thisMinStr = ['0',thisMinStr]; end
    thisSecStr = sprintf('%.2f',thisTimeSec);  while length(thisSecStr) < 5, thisSecStr = ['0',thisSecStr]; end
    timeStr = sprintf('%s:%s:%s',thisHrStr,thisMinStr,thisSecStr);
end

% The below function grabs the sum of the peirwise distances amongst
% members of a cluster
function newD = dFromDist(dist, indices)
    % Nice, fast version by making a submatrix...
    newMat = dist(indices,indices);
    newD = nansum(newMat(:));
end

function w = getW(clusts,dist)
    uClusts = unique(clusts(~isnan(clusts)));
    allD = nan(1,length(uClusts));
    allN = nan(1,length(uClusts));
    for ic = 1:length(uClusts)
        clustInds = find(clusts==uClusts(ic));
        allD(ic) = dFromDist(dist,clustInds);
        allN(ic) = length(clustInds);
    end
    w = sum(allD./(2.*allN));
end


end