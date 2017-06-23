%% This, like klLPP, aims to make Chen & Cai's LSC more user friendly and include the gap

function [outK, obsLabels, gap, gapErr] = klLSCwGap(inMat,varargin)

% Set defaults
thisK=1:6;
randReps = 10;
randType = 'gauss';
options = struct;
gapPerc = 0;

% Decode varargin
varStrInd=find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-p','p'},
            options.p = varargin{varStrInd(iv)+1};
        case {'perc'},
            gapPerc = varargin{varStrInd(iv)+1};
    end
end

%% Let's get started: make a distance matrix
distMat = EuDist2(inMat);

%% Now run the LSC
fprintf('\tDoing Actual Observations...\n');

[obsLabels,obsDists] = doLSC(inMat,thisK);
    

%% Now we generate nReps sets of random data and repeat above analysis

for irand = 1:randReps,
    if mod(irand,5)==1,
        fprintf('\tDoing Random Obs %d of %d...\n',irand,randReps);
    end
    
    switch randType
        case 'unif'
            randVals = min(inMat,[],1)+(max(inMat,[],1)-min(inMat,[],1));
            [~,randDists(irand,:)] = doLSC(randVals,thisK);
        case 'gauss'
            randVals = randn(size(inMat,1),size(inMat,2)).*repmat(nanstd(inMat,[],1),size(inMat,1),1)+repmat(nanmean(inMat,1),size(inMat,1),1);
            [~,randDists(irand,:)] = doLSC(randVals,thisK);
    end
    
end

% randNorm = log(randDists)-repmat(log(randDists(:,1)),1,size(randDists,2));
% obsNorm = log(obsDists)-log(obsDists(1));

randNorm = (log(randDists))-(repmat(log(randDists(:,1)),1,size(randDists,2))-repmat(log(obsDists(1)),size(randDists,1),size(randDists,2)));
obsNorm = log(obsDists);

if gapPerc,
    randNorm = randNorm./repmat(randNorm(:,1),1,size(randNorm,2));
    obsNorm = obsNorm./obsNorm(1);
end

wStar = nanmean(randNorm,1);
sdk = nanstd(randNorm,[],1); %nanmean((randNorm-repmat(wStar,size(randNorm,1),1)).^2);
sk = sdk.*sqrt(1+(1/randReps));

gap = wStar-obsNorm;
gapErr = sk;
gapComp = [gap(2:end)-gapErr(2:end),nan];

outK = find(gap > gapComp,1);
if isempty(outK), outK = length(gap); end;

    function [labels,withinClustDist] = doLSC(inData,k)
        withinClustDist = zeros(1,length(k));
        for ik = 1:length(k),
            labels(:,ik) = LSC(inData,k(ik),options);

            % Loop through each label to get within-cluster distance
            for iik = 1:k(ik),
                centroid = nanmean(inData(labels(:,ik)==iik,:),1);
                thisClustDistMat = EuDist2(inData(labels(:,ik)==iik,:),centroid);
                withinClustDist(ik) = withinClustDist(ik)+nansum(thisClustDistMat);

            end
        end
    end

    function [outIDs,sumDists] = getSortDist(input,myK,kDist)

        nReps = 50;
        kDists = nan(6);
        for ii=myK,
            clear idx sumd
        %     fprintf('Doing k=%d\n',k);
            sumd = nan(ii,nReps);
            for iii = 1:nReps,
                [idx(:,iii),~,sumd(:,iii)] = kmeans(input,ii,'Emptyaction','drop','Distance',kDist);
            end
            [outIDs(:,ii),modDists] = klModeKmeans(idx,sumd);
            [sumDists(ii)] = nansum(modDists);
        end

    end
end