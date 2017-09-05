function [center, dist] = klCluster(vals)
global nClusts stepAlph nit colors seed reps plotClusts

% Set defaults
nClusts = 64;
stepAlph = .05;
nit = 500;
reps = 1;
colors = 'rgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmkyrgbcmky';
seed = 'random';
plotClusts = 1;

% Rotate vals if necessary
if size(vals,1) ~= length(vals),
    vals = vals';
end

% Get the starting variance and centers
rawVar = var(vals,1);
rawCent = nanmean(vals,1);


% Loop through replicates to test for robustness of clusters
for ir = 1:reps,
    [clust(:,ir),cent{ir},dist{ir},mnDist{ir}] = getClust(vals);
end

% Test pairwise whether clusters run together
clustP = zeros(size(clust,1),size(clust,1));
for im = 1:size(clust,1)
    for in = im:size(clust,1),
        clustP(im,in) = sum(clust(im,:) == clust(in,:))/size(clust,2);
    end
end
pTri = tril(clustP);
clustP(pTri ~= 0) = pTri(pTri ~= 0);

sortVect = symrcm(clustP);
sortP = clustP(sortVect,sortVect);

keyboard

    function [outClust, outCent, outDist, outMnDist] = getClust(vals)

    % Get data range and partitions
    sortVals = sort(vals,1);
    minVals = sortVals(1,:);
    for iv = 1:size(vals,2),
        maxVals = sortVals(find(~isnan(sortVals(:,iv)),1,'last'),iv);
    end

    % Get starting point for boundaries of clusters
    switch seed
        case 'random',
            startCents = rand(nClusts,size(vals,2)).*repmat(maxVals-minVals,nClusts,1)+repmat(minVals,nClusts,1);
        otherwise
            startBounds = sortVals(1,:);
            for id = 1:size(vals,2)
                for ic = 1:nClusts,
                    startBounds(ic+1,id) = sortVals(ceil(ic*size(vals,1)*id/nClusts),id);
                    startCents(ic,id) = (startBounds(ic+1,id)+startBounds(ic,id))./2;

                end
            %     keyboard
            end

    end
    % Loop through the clusters to calculate each point's distance to the
    % respective centers to find a closest center

    for ic = 1:nClusts
        rawDiff{ic} = vals - repmat(startCents(ic,:),size(vals,1),1);
        thisDist(:,ic) = sqrt(nansum(rawDiff{ic}.^2,2));
    end

    % Find minimum cluster distance
    for iv = 1:size(vals,1),
        cInd(iv) = find(thisDist(iv,:) == min(thisDist(iv,:)),1);
    end

    % Get mean distances to cluster centers
    for ic = 1:nClusts,
        mnDistRaw(ic) = nanmean(thisDist(cInd == ic,ic),1);
    end

    % Do an initial shift: 10% of the within cluster range?
    for ic = 1:nClusts
        switch seed
            case 'random'
                shiftCents(ic,:) = rand(1,size(vals,2)).*(maxVals-minVals)+minVals;
            otherwise
                shiftCents(ic,:) = startCents(ic,:) - (startBounds(ic+1,:)-startBounds(ic,:)).*.1;
        end
        shiftDiff{ic} = vals-repmat(shiftCents(ic,:),size(vals,1),1);
        shiftDist(:,ic) = sqrt(nansum(shiftDiff{ic}.^2,2));    
    end

    for iv = 1:size(vals,1),
        cShift(iv) = find(shiftDist(iv,:) == min(shiftDist(iv,:)),1);
    end
    distDiff = cellfun(@minus,shiftDiff,rawDiff,'UniformOutput',0);

    for ic = 1:nClusts,
        mnDistShift(ic) = nanmean(shiftDist(cShift == ic,ic),1);
        mnDiff(ic,:) = nanmean(distDiff{ic}(cShift == ic,:),1);
    end

    % Now use the difference in mean distance to get the shift amount for the
    % new clusters. Then iterate over a large number of iterations to achieve
    % convergence

    oldCent = shiftCents;
    oldMnDist = mnDistShift;
    oldDist = shiftDist;
    oMnDiff = mnDistShift;
    oldDiff = shiftDiff;
    for ii = 1:nit
        shiftAmnt = stepAlph.*mnDiff;
        newCent = oldCent+shiftAmnt;
        for ic = 1:nClusts,
            newDiff{ic} = vals-repmat(newCent(ic,:),size(vals,1),1);
            newDist(:,ic) = sqrt(nansum(newDiff{ic}.^2,2));
        end
        for iv = 1:size(vals,1),
            whichClust(iv) = find(newDist(iv,:) ==  min(newDist(iv,:)),1);
        end
        distDiff = cellfun(@minus,newDiff,oldDiff,'UniformOutput',0);
        for ic = 1:nClusts,
            %newMnDist(ic) = nanmean(newDist(whichClust == ic,ic),1);
            newMnDist(ic) = nanmean(newDist(whichClust == ic,ic),1);
    %         mnDiff(ic) = nanmean(distDiff(whichClust == ic,ic),1);
            mnDiff(ic,:) = nanmean(newDiff{ic}(whichClust == ic,:),1);
            nSumDist(ic,:) = nansum(distDiff{ic}(whichClust == ic,:),1);
        end
        distDiff = newMnDist-oldMnDist;
        trackDiff(ii,:) = distDiff;
        oldCent = newCent;
        oldMnDist = newMnDist;
        oldDist = newDist;
    %     for ic = 1:nClusts,
    %         scatter(vals(whichClust == ic),ones(sum(whichClust == ic),1),colors(ic));
    %         vline(oldCent(ic));
    %         hold on;
    %     end
    %     
    %     keyboard
    end

    outCent = oldCent;
    outDist = distDiff;
    outClust = whichClust;
    outMnDist = newMnDist;

    if plotClusts
        figure();
        for ic = 1:nClusts
            if size(vals,2) == 1,
                scatter(vals(outClust == ic),ones(sum(outClust == ic),1),colors(ic));
                hold on;
                vline(outCent(ic));
            elseif size(vals,2) == 2,
                scatter(vals(outClust == ic,1),vals(outClust == ic,2),colors(ic),'filled');
                hold on;
                scatter(outCent(ic,1),outCent(ic,2),'k','filled');

                xErr = nanstd(vals(outClust == ic,1));
                yErr = nanstd(vals(outClust == ic,2));
                errEllipse = klMakeEllipse(xErr,yErr,'-c',[outCent(ic,1),outCent(ic,2)]);
                errSurf = patch(errEllipse(:,1),errEllipse(:,2),'k');
                set(errSurf,'facealpha',.2);
            elseif size(vals,2) == 3,
                scatter3(vals(outClust==ic,1),vals(outClust==ic,2),vals(outClust==ic,3),colors(ic));
                hold on;
                scatter3(outCent(ic,1),outCent(ic,2),outCent(ic,3),'k','filled');

            end
        end
    end
    end
end