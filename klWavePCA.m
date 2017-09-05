monk = {'Gauss','Helmholtz'};
xlFile = 'klDataBookKeeping_mg.xlsx';
types  = {'vis','mov','vismov','none'};

wvSampFreq = 40000;
wvTimeInc  = 1/wvSampFreq;
wvTimes = (1:32).*wvTimeInc.*(10^6);

close all
global excelNum excelAll

% Load column map
makeColLookup
load xlCols

% Initialize some variables
allWaves = []; monkID = []; catFacts = []; catCrit = [];
for im = 1:length(monk),
    [excelNum,~,excelAll] = xlsread(xlFile,monk{im});
    
    myCrit      = strcmp(excelAll(:,col.area),'FEF');
    critInds    = find(myCrit);
    
    % Get waveforms
    monkWaves = klGetWaveforms(critInds,'-m',monk{im});
    allWaves  = cat(1,allWaves,monkWaves);
    
    % Get column map for spike measures
    wdCol       = col.spkStart + (find(strcmpi(col.spkNames,'spkWidth'),1)-1);
    tfrCol      = col.spkStart + (find(strcmpi(col.spkNames,'TfR'),1)-1);
    ampCol      = col.spkStart + (find(strcmpi(col.spkNames,'amplitude'),1)-1);
    relAmpCol   = ampCol + 1;
    mnRateCol   = col.spkStart + (find(strcmpi(col.spkNames,'meanRate'),1)-1);
    fanoCol     = col.spkStart + (find(strcmpi(col.spkNames,'fano'),1)-1);
    cvCol       = col.spkStart + (find(strcmpi(col.spkNames,'cv'),1)-1);
    cv2Col      = col.spkStart + (find(strcmpi(col.spkNames,'cv2'),1)-1);
    lvCol       = col.spkStart + (find(strcmpi(col.spkNames,'lv'),1)-1);
    lvrCol      = col.spkStart + (find(strcmpi(col.spkNames,'lvr'),1)-1);
    mnISICol    = col.spkStart + (find(strcmpi(col.spkNames,'meanISI'),1)-1);

    % Get spike measures from excel file
    rawWd       = excelNum(myCrit,wdCol);
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
    
    factMat = [rawPos,rawWd,rawAmp,rawRange,rawRate,rawFano,rawCV,rawCV2,rawLV,rawISI];
    factNames = {'Pos','Width','Amp','Range','MnRate','Fano','CV','CV2','LV','ISI'};
    
    monkID                  = cat(1,monkID,ones(size(factMat,1),1).*im);
    catFacts                = cat(1,catFacts,factMat);
    catCrit                 = cat(1,catCrit,critInds);
end

keyboard