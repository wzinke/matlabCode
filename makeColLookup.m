% Make column headers for new (15-07-21) excel file
function makeColLookup(recSys)

monk = {'Gauss','Helmholtz','Darwin'};
switch recSys,
    case {1,'plex'},
        xlFile = 'klDataBookKeeping_mg.xlsx';
    case {2,'tdt'},
        xlFile = 'klTDTBookKeeping.xlsx';
end

loadMonk = 1;
[excelNum,excelText,excelAll] = xlsread(xlFile,monk{loadMonk});
while size(excelAll,1) < 5,
    loadMonk = loadMonk+1;
    [excelNum,excelText,excelAll] = xlsread(xlFile,monk{loadMonk});
end

hRow = 4;%find(strcmpi(excelAll(:,2),'Name'),1);
col.sess        = find(strcmpi(excelAll(hRow,:),{'Name'}) | strcmpi(excelAll(hRow,:),{'Session'}),1);
col.count       = find(strcmpi(excelAll(hRow,:),{'Count'}),1);
col.file        = find(strcmpi(excelAll(hRow,:),{'Session'}) | strcmpi(excelAll(hRow,:),{'MG FileName'}),1);
col.chan        = find(strcmpi(excelAll(hRow,:),'chanCode'),1);
col.type        = find(strcmpi(excelAll(hRow,:),'type'),1);
col.typeAlt     = find(strcmpi(excelAll(hRow,:),'typeAlt'),1);
col.area        = find(strcmpi(excelAll(hRow,:),'area'),1);
col.depthChan   = find(strcmpi(excelAll(hRow,:),'depth (chan#)'),1);
col.depth       = find(strcmpi(excelAll(hRow,:),{'depth (um)'}) | strcmpi(excelAll(hRow,:),{'Depth'}),1);
col.spkShape    = find(strcmpi(excelAll(hRow,:),{'wvWidth'}) | strcmpi(excelAll(hRow,:),{'spkWidth'}),1);
col.shapeNames  = excelAll(hRow,col.spkShape:find(strcmpi(excelAll(hRow,:),'normRange')));
spkStarts       = find(strcmpi(excelAll(hRow,:),'meanRate'),2);
spkEnds         = find(strcmpi(excelAll(hRow,:),' < 3ms'),2);
col.spkStartT   = spkStarts(1);
col.spkNamesT   = excelAll(hRow,spkStarts(1):spkEnds(1));
col.spkStartB   = spkStarts(2);
col.spkNamesB   = excelAll(hRow,col.spkStartB:spkEnds(2));
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
col.SNR 		= find(strcmpi(excelAll(hRow,:),'SNR'),1);

switch recSys,
    case {1,'plex'},
        plexCols = col;
        save plexCols.mat plexCols
    case {2,'tdt'},
        tdtCols = col;
        save tdtCols.mat tdtCols
end
