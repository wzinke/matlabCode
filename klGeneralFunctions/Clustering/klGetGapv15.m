%% This function aims to recreate the "gap" statistic for cluster quality
%  Tibshirani et al 2001 describe gap statistic for identifying number of
%  clusters in a given data set

%  v8 removes 'boot' randType ability for ease of K-Means distances for
%  speed

function [kHat, gapVect, s, summaryStruct] = klGetGapv10(obs,groups,varargin)

% Set defaults
type = 'euc';
isSim = 0;
catVars = zeros(1,size(obs,2));
simMat = [];
doN = 1;
randSize = 100;
subSampPerc = .1;
randReps = 40;
randType = 'rand';
getPCA = 0;
reClustType = 'agglom';
randValType = 'pca';
pcaRandType = 'gauss';
subSetPerc = .5;
minConsec = 1;
doPerc = 0;
nS = 1;
dim = 1;
weights = ones(1,size(obs,2));

% Set reclustering defaults
nMembers = 10;
nMembersRaw = nMembers;
refine = 1;
refType = 'same';
kType = 'correlation';
kReps = 30;
goodType = 'sumsquares';
nPrint = 100;
print = 0;
replace = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','type'}
            type = varargin{varStrInd(iv)+1};
        case {'-s','sim'}
            simMat = varargin{varStrInd(iv)+1};
            isSim = 1;
        case {'cat'}
            catVars = varargin{varStrInd(iv)+1};
        case {'-n','n'},
            nMembers = varargin{varStrInd(iv)+1};
        case {'-r','randtype'},
            randType = varargin{varStrInd(iv)+1};
        case {'reclust'},
            reClustType = varargin{varStrInd(iv)+1};
        case {'randreps'},
            randReps = varargin{varStrInd(iv)+1};
        case {'randval'},
            randValType = varargin{varStrInd(iv)+1};
        case {'replace'},
            replace = varargin{varStrInd(iv)+1};
        case {'-c','consec','minconsec'},
            minConsec = varargin{varStrInd(iv)+1};
        case {'ktype'},
            kType = varargin{varStrInd(iv)+1};
        case {'-d'},
            cDists = varargin{varStrInd(iv)+1};
        case {'perc'},
            doPerc = varargin{varStrInd(iv)+1};
        case {'nsig'},
            nS = varargin{varStrInd(iv)+1};
        case {'dim'},
            dim = varargin{varStrInd(iv)+1};
        case {'-w','weight'},
            weights = varargin{varStrInd(iv)+1};
        case {'refine'},
            refine = varargin{varStrInd(iv)+1};
    end
end

obsRaw = obs;
catVarInds = find(catVars);
nonCatInds = find(~catVars);

if ismember(type,{'prod'}),
    maxObs = nanmax(obs,[],1);
    minObs = nanmin(obs,[],1);
    normObs = (obs-repmat(minObs,size(obs,1),1))./repmat((maxObs-minObs),size(obs,1),1);
    normRaw = normObs;
    obs = normObs;
end

% Make similarity matrix, if necessary
% Start by computing pair-wise distances
if ~exist('cDists','var')
    fprintf('Initializing pairwise distances...')
    if ismember(type,{'corr'}),
        % Clean up obs, if "corr"
        goodCols = zeros(1,size(obs,2));
        for ic = 1:size(obs,2),
            goodCols(ic) = sum(~isnan(obs(:,ic))) == size(obs,1);
        end
        obs = obs(:,logical(goodCols));
        catVars = catVars(logical(goodCols));
        simMat = 1-corr(obs',obs');
    else
        simMat = nan(size(obs,1),size(obs,1));
        for im = 1:size(obs,1),
            for in = (im+1):size(obs,1),
                simMat(im,in) = klGetDistv2(obs(im,:),obs(in,:),'-t',type,'cat',catVars,'-w',weights);
%                 if simMat(im,in) > 10000,
%                     keyboard
%                 end
                simMat(in,im) = simMat(im,in);
            end
        end
    end
    fprintf('Done!\n');
end

% if ismember(randType,{'rand'}),
%     randObs = nan(size(obs,1),size(obs,2),randSize);
%     for iReps = 1:ceil(randSize/subSampPerc),
%         for iCol = 1:size(obs,2),
%             randObs(:,iCol,iReps) = (rand(size(obs,1),1).*range(obs(:,iCol)))-min(obs(:,iCol));
%             if catVars(iCol),
%                 randObs(:,iCol,iReps) = round(randObs(:,iCol,iReps));
%             end
%         end
%         for ir = 1:size(obs,1),
%             randSimMat(ir,ir,iReps) = 0;
%             for ic = ir:size(obs,1),
%                 thisDist = klGetDist(randObs(ir,:,iReps),randObs(ic,:,iReps),'-t',type,'cat',catVars);
%                 randSimMat(ir,ic,iReps) = thisDist;
%                 randSimMat(ic,iReps) = thisDist;
%             end
%         end
%     end
% end

if getPCA || ismember(randValType,{'pca'}),
    [coeffs, score] = pca(obs(:,~catVars));%,'Centered',0);
end

% Start k loop
for ik = 1:size(groups,2),
    clusts = unique(groups(:,ik)); clusts(isnan(clusts)) = [];
    k(ik) = max(clusts);
%     fprintf('\tGetting gap for k=%d\n',k(ik));
    
    % Definitions:
    % Dr = sum(d(ij)) where i, j are in cluster r
    % Wk = sum((1/2nr)*Dr)
    
    if ~exist('cDists','var'),
        D = zeros(1,length(clusts));
        meanRandD = zeros(1,length(clusts));
        n = zeros(1,length(clusts));
        randD = zeros(randSize,length(clusts));
        shuffGroups = groups(randperm(1:size(groups,1)),length(clusts));
        
        % Start cluster loop
        for ic = 1:max(clusts),
            clustInds = find(groups(:,ik) == clusts(ic));
            shuffInds = find(shuffGroups(:,ik)==clusts(ik));
            
            n(ic) = length(clustInds);
            % Start within cluster loop
            D(ic) = dFromSim(simMat,clustInds);

            % Here is where I differ (for now):
            % As a reference I intend to use a randomization procedure (in the
            % future I will rearrange such that it can do this or use the
            % uniform distribution described in Tibshirani et al 2001
            switch randType,
                case {'subsample','bootstrap','boot'},
                    for irand = 1:ceil(randSize/subSampPerc),
                        if replace,
                            randInds = randi([1,size(groups,1)],[1,size(groups,1)]);
                        else
                            randInds = randperm(size(groups,1));
                        end
                        randD(irand,ic) = dFromSim(simMat,randInds(1:length(clustInds)));
                    end
                case {'rand'},
    %                 for irand = 1:ceil(randSize/subSampPerc),
    %                     randD(irand,ic) = dFromSim(randSimMat,clustInds);
    %                 end
                case {'pca'}
                    for iCol = 1:size(obs,2),
                    end
                case {'shuffle'},
%                     randD(iarand,ic) = dFromSim
            end

        end
    else
        D = cDists(:,ik); D(isnan(D)) = [];
        uGrps = unique(groups(:,ik));
        for ig = 1:length(uGrps),
            n(ig,1) = sum(groups(:,ik) == uGrps(ig));
        end
    end
    W(ik) = sum(D./(2.*n));
    sumD(ik) = sum(D);
    
end

if ismember(randType,{'rand'}),
    fprintf('\tGetting %d random observations...',randReps);
    for irand = 1:randReps,
        fprintf('%d',irand);
%         fprintf('Getting random observation %d of %d\n',irand,randReps);
        % Get random observations
        switch randValType,
            case {'uniform','unif'},
                % Get uniform values
                unifMat = rand(size(obs));
                
                % Get range, min, max
                minMat = repmat(min(obs,[],1),size(obs,1),1);
                maxMat = repmat(max(obs,[],1),size(obs,1),1);
                rangeMat = maxMat-minMat;
                
                % Get random observations
                randObs = (unifMat.*rangeMat)+minMat;
            case {'gauss'},
                % Get uniform values
                unifMat = randn(size(obs));
                % Get random observations
                randObs = (unifMat.*repmat(nanstd(obs,[],1),size(obs,1),1))+repmat(nanmean(obs,1),size(obs,1),1); 
            case {'pca'},
                randObs = nan(size(obs));
                
                % Get uniform values
                 
                % Get range, min, max
                switch pcaRandType
                    case 'unif',
                        unifMat = rand(size(score));
                        minMat = repmat(min(score,[],1),size(score,1),1);
                        maxMat = repmat(max(score,[],1),size(score,1),1);
                        rangeMat = maxMat-minMat;
                
                        % Get random observations
                        randScore = (unifMat.*rangeMat)+minMat;
                        randObs(:,~catVars) = (randScore*coeffs');
                    case 'gauss',
                        unifMat = randn(size(score));
                        randObs = (unifMat.*repmat(nanstd(score,[],1),size(score,1),1))+repmat(nanmean(score,1),size(score,1),1); 
                end
                % Replace column of randObs with uniformally picked
                % categorical option from original observations
                for iv = 1:length(catVarInds),
                    uCat = unique(obs(:,catVarInds(iv)));
                    randObs(:,catVarInds(iv)) = uCat(randi([1,length(uCat)],[size(obs,1),1]));
                end
                
            case {'subset'},
                randObs = obs(randi([1,size(obs,1)],[ceil(size(obs,1).*subSetPerc),1]),:);
                nMembers = ceil(nMembersRaw*subSetPerc);
        end
        
        switch reClustType
            case 'agglom'
                [randClustIDs,~,randDistMat] = klAgglomv10(randObs,'-k',k,'-r',refine,'-p',print,'-t',type,'-n',nMembers,'refType',refType,'kType',kType,'nprint',nPrint);
                randMat = zeros(size(randDistMat));
                for ir = 1:size(randDistMat,1),
                    for ic = (ir+1):size(randDistMat,2),
                        randMat(ir,ic) = randDistMat(ir,ic);
                        randMat(ic,ir) = randDistMat(ir,ic);
                    end
                end
            case 'kmeans'
                warning off 
                for ir = 1:size(randObs,1),
                    firstVal(ir) = find(~isnan(randObs(ir,:)),1);
                    lastVal(ir) = find(~isnan(randObs(ir,:)),1,'last');
                end
                clpRand = randObs(:,max(firstVal):min(lastVal));
                randKDist = nan(length(k));
                for ik = 1:length(k),
                    nTries = 0;
                    randClustIDs(:,ik) = nan(size(obs,1),1);
                    while sum(isnan(randClustIDs(:,ik))) && nTries < 25,
                        try
                            [randClustIDs(:,ik),~,randKDist(1:ik,ik)] = kmeans(clpRand,k(ik),'Distance',kType,'Replicates',kReps,'EmptyAction','drop');
                        catch
                            randClustIDs(:,ik) = nan(size(obs,1),1);
                            randKDist(1:ik,ik) = nan(1,ik);
                        end
                        nTries = nTries+1;
                    end
                    if length(unique(randClustIDs(:,ik))) > k(ik),
                        keyboard
                    end
                end
%                 if ismember(type,{'corr'}),
%                     % Clean up obs, if "corr"
%                     goodCols = zeros(1,size(clpRand,2));
%                     for ic = 1:size(clpRand,2),
%                         goodCols(ic) = sum(~isnan(clpRand(:,ic))) == size(obs,1);
%                     end
%                     clpRand = clpRand(:,logical(goodCols));
%                     randMat = 1-corr(clpRand',clpRand');
%                 else
%                     randMat = nan(size(clpRand,1),size(clpRand,1));
%                     for im = 1:size(clpRand,1),
%                         for in = (im+1):size(clpRand,1),
%                             randMat(im,in) = klGetDist(clpRand(im,:),clpRand(in,:),'-t',type,'cat',catVars);
%                             randMat(in,im) = randMat(im,in);
%                         end
%                     end
%                 end

            case 'pca'
                warning off
                for ir = 1:size(allWaves,1),
                    firstVal(ir) = find(~isnan(randObs(ir,:)),1);
                    lastVal(ir) = find(~isnan(randObs(ir,:)),1,'last');
                end
                clpRand = randObs(:,max(firstVal):min(lastVal));
                [~,scores,~,~,expl] = pca(clpRand);
                cumExpl = cumsum(expl);
                nDims = find(cumExpl > 95,1);
                switch pcaTypeWv,
                    case 'kmeans',
                        for ik = 1:length(k),
                            for ir = 1:nReps,
                                [pcaIDX(:,ir)] = kmeans(scores(:,1:nDims),k(ik));
                            end
                            randClustIDs(:,ik) = klModeKmeans(pcaIDX);
                        end
                        simMat = [];
                    case 'agglom'
                        randClustIDs = klAgglomv10(scores,'-k',k,'-r',refine,'-p',print,'-t',pcaAgglomType,'-n',nWvMembers,'refType',refType,'kType',kTypeWv,'nprint',nPrint,'-w',weights);
                    otherwise
                        randClustIDs = klAgglomv10(scores,'-k',k,'-r',refine,'-p',print,'-t',pcaAgglomType,'-n',nWvMembers,'refType',pcaRefType,'kType',kTypeWv,'nprint',nPrint,'-w',weights);
                end
            otherwise
                [randClustIDs,~,randDistMat] = klAgglomv10(randObs,'-k',k,'-r',refine,'-p',print,'-t',type,'-n',nMembers,'refType',refType,'kType',kType,'nprint',nPrint);
                randMat = zeros(size(randDistMat));
                for ir = 1:size(randDistMat,1),
                    for ic = (ir+1):size(randDistMat,2),
                        randMat(ir,ic) = randDistMat(ir,ic);
                        randMat(ic,ir) = randDistMat(ir,ic);
                    end
                end
        end

        % Now, get Wstar for this repetition
        for ik = 1:length(k),
%             randClusts = unique(randClustIDs(:,ik)); randClusts(isnan(randClusts)) = [];
            randClusts = 1:ik;
            randN=zeros(1,length(randClusts));
            Dstar = nan(1,length(randClusts));
            for ic = 1:max(randClusts),
                randInds = find(randClustIDs(:,ik) == randClusts(ic));
                randN(ic) = length(randInds);
                if ~strcmp(reClustType,'kmeans'),
                    Dstar(ic) = dFromSim(randMat,randInds);
                end
            end
            if strcmp(reClustType,'kmeans'),
                Dstar(:) = randKDist(1:length(randClusts),ik);
            end
            Wstar(irand,ik) = nansum(Dstar./(2.*randN));
            sumDstar(irand,ik) = nansum(Dstar);
    %             gapVect(irand,ik) = log(Wstar(irand,ik)-log(W(ik)));
            if Wstar(irand,ik) == -inf,
                keyboard
            end
        end
        for ib = 1:length(sprintf('%d',irand)),
            fprintf('\b',irand);
        end
    end
    fprintf('%d\n',irand);
    
    % Does it make sense to normalize W and Wstar by their total starting
    % distances? If "gap" is the "boost" in reduced distance, if the
    % uniformally distributed dataset has more distance to "burn", perhaps
    % it should be thought of as a percentage?
    
    gapVect = nanmean(log(Wstar)-repmat(log(W),size(Wstar,1),1),1);
    gapVect2 = nanmean(log(sumDstar)-repmat(log(sumD),size(sumDstar,1),1),1);
    
    percVect = (log(sumD(1))-log(sumD))./log(sumD(1));
    percVectStar = (log(sumDstar(1))-log(sumDstar))./log(sumDstar(1));
    
    percGap = nanmean(repmat(percVect,size(percVectStar,1),1)-percVectStar,1);
    
    lBar = nanmean(log(Wstar),1);
    lBar2 = nanmean(log(sumDstar),1);
    
    sd = sqrt(nanmean((log(Wstar)-repmat(lBar,size(Wstar,1),1)).^2,1));
    sd2 = sqrt(nanmean((log(sumDstar)-repmat(lBar2,size(sumDstar,1),1)).^2,1));
    s = sd.*sqrt(1+(1/randReps));
    s2 = sd2.*sqrt(1+(1/randReps));
    
    percBar = nanmean(log(percVectStar),1);
    sdPerc = sqrt(nanmean((log(percVectStar)-repmat(percBar,size(percVectStar,1),1)).^2,1));
    sPerc = sdPerc.*sqrt(1+(1/randReps));
    
    
%     gapVect = nanmean(gapVect,1);
%     lBar = nanmean(log(Wstar),1);
%     sd = sqrt(nanmean((log(Wstar)-repmat(lBar,size(Wstar,1),1)).^2,1));
%     s = sd.*sqrt(1+(1/randReps));
    gapCompVect = [gapVect(2:end) - nS.*s(2:end),nan];
    gapComp2 = [gapVect2(2:end) - nS.*s2(2:end),nan];
    percCompVect = [percGap(2:end) - nS.*sPerc(2:end),nan];

    switch dim,
        case {1,'forward','for'},
            findDim = 'first';
            
            critVect = gapVect > gapCompVect;
            consecVect = klGetConsecutive(critVect);

            critVect2 = gapVect > gapComp2;
            consecVect2 = klGetConsecutive(critVect2);

            percCrit = percGap > percCompVect;
        case {2,'backward','back'},
            findDim = 'last';
            critVect = gapVect < gapCompVect;
            consecVect = klGetConsecutive(critVect);

            critVect2 = gapVect < gapComp2;
            consecVect2 = klGetConsecutive(critVect2);

            percCrit = percGap < percCompVect;
        otherwise
            findDim = 'first';
            
            critVect = gapVect > gapCompVect;
            consecVect = klGetConsecutive(critVect);

            critVect2 = gapVect > gapComp2;
            consecVect2 = klGetConsecutive(critVect2);

            percCrit = percGap > percCompVect;
    end
    
    if doPerc,
        kVect = k(find(klGetConsecutive(percCrit) >= minConsec,1,findDim));
    else
        if strcmp(reClustType,'kmeans'),
            kVect = k(find(klGetConsecutive(critVect2) >= minConsec,1,findDim));
        else
            kVect = k(find(klGetConsecutive(critVect) >= minConsec,1,findDim)); 
        end
        if isempty(kVect),
            if dim==2,
                kVect = min(k)-1;
            else
                kVect = max(k);
            end
        end
    end
    if dim == 2,
        kVect = kVect+1;
    end
end

if ~ismember(randType,{'rand'}),
    for irand = 1:randReps,
        % khat = smallest k such that gap(k) >= (gap(k+1)-s(k+1))
        gapCompVect(irand,:) = [gapVect(irand,2:end) - nS.*s(irand,2:end),nan];
        kVect(irand) = k(find(gapVect(irand,:) > gapCompVect(irand,:),1));
        fprintf('*** Update this section of the code with forward/backward stuff...\n');
        keyboard
    end
end

kHist = hist(kVect,k);
kHat = k(kHist == max(kHist));

summaryStruct.kHat = kHat;
summaryStruct.gapVect = gapVect;
summaryStruct.s = s;
summaryStruct.W = W;
summaryStruct.Wstar = Wstar;
summaryStruct.kVect = kVect;

function bigD = dFromSim(similarity, indices)
    bigD = 0;
    for ii = 1:length(indices),
        for iii = (ii+1):length(indices),
            bigD = bigD + similarity(indices(ii),indices(iii));
        end
    end
end

end