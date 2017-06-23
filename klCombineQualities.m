monk = {'Gauss','Helmholtz','Darwin'};
sav = 1;
writ = 1;
xlFile = './klDataBookKeeping_mg.xlsx';

makeColLookup
load xlCols
for im = 1:length(monk),
    clear monkSorts startPt endPt monkQuals
    % Get this monk's sort files
    mySorts = dir(sprintf('./%sQualities*.mat',monk{im}(1)));
    mySortNames = {mySorts.name};
    
    % Loop through each file name and get the end point
    for is = 1:length(mySortNames),
        undInd  = find(ismember(mySortNames{is},'_'),1,'last');
        hyphInd = find(ismember(mySortNames{is},'-'),1,'last');
        suffInd = find(ismember(mySortNames{is},'.'),1,'last');
        
        startPt(is) = str2num(mySortNames{is}((undInd+1):(hyphInd-1)));
        endPt(is) = str2num(mySortNames{is}((hyphInd+1):(suffInd-1)));
    
        load(mySortNames{is});
        allSorts.(monk{im}).sortInfo(startPt(is):endPt(is)) = sortQualities.(monk{im}).sortInfo(startPt(is):endPt(is));
        allSorts.(monk{im}).quals(startPt(is):endPt(is),:) = sortQualities.(monk{im}).quals(startPt(is):endPt(is),:);
    end
    
    if writ,
        xlswrite(xlFile,allSorts.(monk{im}).quals(5:end,:),monk{im},sprintf('%s5',num2abc(col.SNR)));
    end
%     allSorts.(monk{im}).sortInfo = monkSorts;
%     allSorts.(monk{im}).quals = monkQuals;
end

if sav,
    save('allSorts.mat','allSorts');
end