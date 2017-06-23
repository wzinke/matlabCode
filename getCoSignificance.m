monk = {'Gauss','Helmholtz'};

for im = 1:2,

[excelNum,excelText,excelAll] = xlsread(xlFile,monk{im});
excelNum = cat(1,nan(1,size(excelNum,2)),excelNum);
excelNum = cat(2,nan(size(excelNum,1),size(excelAll,2)-size(excelNum,2)),excelNum);
posCol          = find(strcmp(excelAll(1,:),'posMax'));
areaCol         = find(strcmp(excelAll(1,:),'Area'));
statColStart = find(strcmpi(excelAll(1,:),'vTrans'));
statColEnd = find(strcmpi(excelAll(1,:),'Err'));
statTests = excelAll(1,statColStart:statColEnd);

for ia = 1:length(areaList),
    sigErr = strcmpi(excelAll(:,areaCol),areaList{ia}) & excelNum(:,statColEnd) < .05;
    for it = 1:length(statTests)-2
        for ip = 0:1
            otherSigTests(it,(ip+1),ia) = sum(sigErr & excelNum(:,statColStart+(it-1)) < .05 & excelNum(:,posCol) == ip);
            outOf(it,(ip+1),ia)         = sum(sigErr & excelNum(:,posCol) == ip);
        end
    end
    percSigTests(:,:,ia) = otherSigTests(:,:,ia)./outOf(:,:,ia);
    figure(ia+(10*(im-1)));
    bar(percSigTests(:,:,ia)); set(gca,'XTickLabel',statTests(1:(end-2)));
    title(sprintf('%s - Area %s',monk{im},areaList{ia}));
end

end