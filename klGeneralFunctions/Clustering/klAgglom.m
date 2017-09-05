%% klAgglom takes an mxn matrix obs where rows are observations and
%  columns are data points. Note that in comments, the recursion is such
%  that I don't want to distinguish between individuals and clusters, so I
%  consider individuals as one-member-clusters.

function [clusters, icd, clusMembers, linkMat] = klAgglom(obs,varargin)

% Set defaults
type        = 'corr';
combType    = 'all'; % use "mean" or "all"
k           = 6;
nMembers    = 10;
refine      = 1;
nPrint      = 50;
print       = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','type'},
            type = varargin{varStrInd(iv)+1};
        case {'-k','k'},
            k = varargin{varStrInd(iv)+1};
        case {'-r','refine'},
            refine = varargin{varStrInd(iv)+1};
        case {'-p','print'},
            print = 1;
        case {'-c','comb'},
            combType = varargin{varStrInd(iv)+1};
        case {'-n'}
            nMembers = varargin{varStrInd(iv)+1};
    end
end

% Start by computing pair-wise distances
distMat = nan(size(obs,1),size(obs,1));
for im = 1:size(obs,1),
    for in = (im+1):size(obs,1),
        switch type,
            case 'corr',
                thisDist = 1-corr(obs(im,:)',obs(in,:)');
            case 'euc',
                thisDist = sqrt(nansum((obs(im,:)-obs(in,:)).^2));
        end
        distMat(im,in) = thisDist;
    end
end

% Save distance and observations in raw form in case I manipulate them in
% the next step...
distRaw = distMat;
obsRaw  = obs;

% Initialize clusMembers
clusMembers = num2cell((1:size(distRaw,1))');
linkMat = [];
linkIDs(1:size(distRaw,1),1) = 1:size(distRaw,1);
nLoops = 0;
while size(distMat,1) > 1,
    if mod(size(distMat,1),nPrint) == 0 && print,
        fprintf('%d members remaining...\n',size(distMat,1));
    end
    nLoops = nLoops+1;
        % Get the row and column corresponding to the minimum distance
        % remaining. If there are multiple, take just the first index
        % because other indices should be captured in the next loop
        % through.
        
        [minRow,minCol] = find(distMat == min(distMat(:)),1);
        nRow = length(clusMembers{minRow,end});
        nCol = length(clusMembers{minCol,end});
        
        % Combine observations (weighted mean)
        combObs = sum([obs(minRow,:).*nRow; obs(minCol,:).*nCol],1)./(nRow+nCol);
        
        % Get the indices of clusMembers to keep
        nonClust = 1:size(distMat,1); nonClust(ismember(nonClust,[minRow,minCol])) = [];
        newMems = clusMembers(nonClust,end);
        newMems{end+1} = [clusMembers{minRow,end},clusMembers{minCol,end}];
        theseMems = newMems{end};
        
        % Let's make a linkage matrix. Remember:
        % For example, suppose there are 30 initial nodes and at 
        % step 12 cluster 5 and cluster 7 are combined. Suppose their 
        % distance at that time is 1.5. Then Z(12,:) will be [5, 7, 1.5]. 
        % The newly formed cluster will have index 12 + 30 = 42. If 
        % cluster 42 appears in a later row, it means the cluster created 
        % at step 12 is being combined into some larger cluster.
        linkMat = cat(1,linkMat,[linkIDs(minRow),linkIDs(minCol),min(distMat(:))]);
%         newLink = linkIDs(minRow)+linkIDs(minCol)+size(distRaw,1);
        newLink = nLoops + size(distRaw,1);
        linkIDs = linkIDs(nonClust);
        linkIDs(end+1) = newLink;
        
        % Debugging to be sure we're not overcounting, should be worked
        % out...
        if any(ismember(clusMembers{minRow,end},clusMembers{minCol,end})),
            keyboard
        end
        
        % Save combined members in clusMembers
        clusMembers(1:length(newMems),size(clusMembers,2)+1) = newMems;
        
        % Put new cluster distances at the end of "obs" after removing the
        % members of that cluster
        obs = obs(nonClust,:);
        obs = cat(1,obs,combObs);
        
        % Cut the cluster members from distMat
        distMat = distMat(nonClust,nonClust);
        
        % Calculate distance between other clusters and the most recent
        % combination
        distMat = cat(1,distMat,nan(1,size(distMat,2)));
        distMat = cat(2,distMat,nan(size(distMat,1),1));
        for im = 1:(size(distMat,1)-1),
            switch combType
                case 'mean'
                    switch type,
                        case 'corr',
                            thisDist = 1-corr(obs(im,:)',obs(end,:)');
                        case 'euc',
                            thisDist = sqrt(nansum((obs(im,:)-obs(end,:)).^2));
                    end
                case 'all'
                    tstMems = clusMembers{im,end};
                    allCorrs = distRaw(tstMems,theseMems);
                    allCorrsTrans = distRaw(theseMems,tstMems)';
                    noNanCorrs = nan(size(allCorrs));
                    noNanCorrs(~isnan(allCorrs)) = allCorrs(~isnan(allCorrs));
                    noNanCorrs(~isnan(allCorrsTrans)) = allCorrsTrans(~isnan(allCorrsTrans));

                    thisDist = nanmean(noNanCorrs(:));
            end
            distMat(end,im) = thisDist;
        end
        
        % Let's make a linkage matrix. Remember:
        % For example, suppose there are 30 initial nodes and at 
        % step 12 cluster 5 and cluster 7 are combined. Suppose their 
        % distance at that time is 1.5. Then Z(12,:) will be [5, 7, 1.5]. 
        % The newly formed cluster will have index 12 + 30 = 42. If 
        % cluster 42 appears in a later row, it means the cluster created 
        % at step 12 is being combined into some larger cluster.
                
end



% Let's put the next part in a loop so we can ask for multiple "k" values
% without having to rerun the time consuming correlation part above.
clusters = nan(size(obsRaw,1),length(k));
        
for ik = 1:length(k),
    % Go backwards and find the last instance with "k" clusters of "nMembers"
    % or more
    nClusts = 1;
    nLoops = 0;
    premExit = 0;

    while nClusts < k(ik)
        nLoops = nLoops + 1;
        clustLens = cellfun(@length,clusMembers(:,end-nLoops));
        nClusts = sum(clustLens >= nMembers);
        if nClusts == 0,
            fprintf('Criteria not reached before waveforms are too scattered. Exiting loop\n');
            premExit = 1;
            continue;
        end
    end

    if premExit,
        clusters = nan(size(obsRaw,1),1);
        return;
    end

    % Get mean data points by cluster
    clustInds   = find(clustLens >= nMembers);
    clustersPre = nan(size(obsRaw,1),1);
    for ic = 1:length(clustInds),
        memInds = clusMembers{clustInds(ic),end-nLoops};
        clustMeans(ic,:) = nanmean(obsRaw(memInds,:),1);
        clustersPre(memInds) = ic;
    end

    % If asked to, loop through all observations and compare to the mean data
    % points and reassign each observation based on the nearest distance
    if refine
        for io = 1:size(obsRaw,1),
            testDist = nan(length(clustInds),1);
            for ic = 1:length(clustInds),
                switch type,
                    case 'corr',
                        testDist(ic) = 1-corr(obsRaw(io,:)',clustMeans(ic,:)');
                    case 'euc',
                        testDist(ic) = sqrt(nansum((obsRaw(io,:)-clustMeans(ic,:)).^2));
                end
            end
            % I haven't figured out how to deal with identical distances yet...
            % put in a catch for debugging this when necessary
            if sum(testDist == min(testDist)) > 1,
                keyboard
            end
            clusters(io,ik) = find(testDist == min(testDist));
            clear testDist;
        end
    else
        clusters(:,ik) = clustersPre;
%         [~,clusters(:,ik)] = dendrogram(linkMat,k(ik),'display','off');
        
    end

    % Get final group means, then calculate inter-group differences
    uClust = unique(clusters(~isnan(clusters(:,ik)),ik));
    for ic = 1:length(uClust),
        clustMn(ic,:) = nanmean(obsRaw(clusters(:,ik) == uClust(ic),:),1);
    end
    icd{ik} = nan(length(uClust),length(uClust));
    for im = 1:length(uClust),
        for in = im:length(uClust),
            switch type,
                case 'corr',
                    thisDist = 1-corr(clustMn(im,:)',clustMn(in,:)');
                case 'euc',
                    thisDist = sqrt(nansum((clustMn(im,:)-clustMn(in,:)).^2));
            end
            icd{ik}(im,in) = thisDist;
        end
    end
end