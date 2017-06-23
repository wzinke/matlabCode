% Make column headers for new (15-07-21) excel file
function makeColLookup()

monk = {'Gauss','Helmholtz'};
xlFile = './klDataBookKeeping_mg.xlsx';

[excelNum,excelText,excelAll] = xlsread(xlFile,monk{1});

hRow = find(strcmpi(excelAll(:,2),'Name'),1);
col.sess        = find(strcmpi(excelAll(hRow,:),'Name'),1);
col.count       = find(strcmpi(excelAll(hRow,:),'Count'),1);
col.file        = find(strcmpi(excelAll(hRow,:),'MG FileName'),1);
col.chan        = find(strcmpi(excelAll(hRow,:),'chanCode'),1);
col.type        = find(strcmpi(excelAll(hRow,:),'type'),1);
col.typeAlt     = find(strcmpi(excelAll(hRow,:),'typeAlt'),1);
col.area        = find(strcmpi(excelAll(hRow,:),'area'),1);
col.depthChan   = find(strcmpi(excelAll(hRow,:),'depth (chan#)'),1);
col.depth       = find(strcmpi(excelAll(hRow,:),'depth (um)'),1);
col.spkStart    = find(strcmpi(excelAll(hRow,:),'spkWidth'),1);
col.spkNames    = excelAll(hRow,col.spkStart:find(strcmpi(excelAll(hRow,:),' < 3ms'),1));
col.mgPvals     = find(strcmpi(excelAll(hRow,:),'vTrans'),1);
col.visIndex    = find(strcmpi(excelAll(hRow,:),'visIndex'),1);
col.movIndex    = find(strcmpi(excelAll(hRow,:),'movIndex'),1);
col.visPoiss    = find(strcmpi(excelAll(hRow,:),'poissVisLat'),1);
col.movPoiss    = find(strcmpi(excelAll(hRow,:),'poissMovLat'),1);
col.visLat      = find(strcmpi(excelAll(hRow,:),'vis2Sig'),1);
col.movLat      = find(strcmpi(excelAll(hRow,:),'mov2Sig'),1);
col.visSupp     = find(strcmpi(excelAll(hRow,:),'visSupp'),1);
col.movSupp     = find(strcmpi(excelAll(hRow,:),'movSupp'),1);
col.vTune       = find(strcmpi(excelAll(hRow,:),'visDir'),1);
col.mTune       = find(strcmpi(excelAll(hRow,:),'movDir'),1);
col.vTST        = find(strcmpi(excelAll(hRow,:),'vLat-Type'),1);
col.mTST        = find(strcmpi(excelAll(hRow,:),'mLat-Type'),1);
save xlCols.mat col