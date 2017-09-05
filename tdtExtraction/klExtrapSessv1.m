function klExtrapSessv1(myFile)

procDir = 'y:/users/kaleb/dataProcessed';

for ic = 1:64,
    fprintf('Extrapolating Channel %d sort...\n\t',ic);
    printStr = 'Loading data...';
    fprintf('%s',printStr);
    load(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',procDir,myFile,ic));
    for ib = 1:length(printStr), fprintf('\b'); end
    printStr = 'Extrapolating Negatives...';
    fprintf('%s',printStr);
    if chanSorts.neg.k > 1,
        chanSorts.neg.allIDX = klExtrapSortv1(chanSorts.neg.allScores(:,1:2),chanSorts.neg.allScores(chanSorts.neg.subInds,1:2),chanSorts.neg.idx(:,chanSorts.neg.k));
    else
        chanSorts.neg.allIDX = uint8(ones(size(chanSorts.neg.allScores,1),1));
    end
    for ib = 1:length(printStr), fprintf('\b'); end
    printStr = 'Extrapolating Positives...';
    fprintf('%s',printStr);
    if chanSorts.pos.k > 1,
        chanSorts.pos.allIDX = klExtrapSortv1(chanSorts.pos.allScores(:,1:2),chanSorts.pos.allScores(chanSorts.pos.subInds,1:2),chanSorts.pos.idx(:,chanSorts.pos.k));
    else
        chanSorts.pos.allIDX = uint8(ones(size(chanSorts.pos.allScores,1),1));
    end
    for ib = 1:length(printStr), fprintf('\b'); end
    printStr = 'Saving...';
    fprintf('%s',printStr);
    save(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',procDir,myFile,ic),'chanSorts','-v7.3');
    for ib = 1:length(printStr), fprintf('\b'); end
    printStr = 'Done!';
    fprintf('%s',printStr);
    fprintf('\n');
end
    