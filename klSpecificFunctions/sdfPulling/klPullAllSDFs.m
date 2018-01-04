function [goodSDF,goodSDFTimes,goodAreas, goodRecSys, goodRow, goodMonks, goodOtherTasks, goodSess, goodChans, goodParams] = klPullAllSDFs(task,varargin)

% Set defaults
areaCrit = {'FEF','SEF','F2','SC'};
doReward = 0;
clip = 0;
recCrit = {'brad','plex','tdt','rich'};

% decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-a','area'}
            areaCrit = varargin{varStrInd(iv)+1};
        case {'-r','reward'}
            doReward = varargin{varStrInd(iv)+1};
        case {'rec'}
            recCrit = varargin{varStrInd(iv)+1};
        case {'-c'}
            clip = varargin{varStrInd(iv)+1};
    end
end

allRecSys = {};

if ismember('brad',recCrit)
    [bradSDF,bradTimes,bradAreas,bradRows,bradMonks,bradOtherTasks,bradSess,bradChans,bradParams] = klPullBradenSEF(task,'-a',areaCrit,'-r',doReward,'-c',clip);
    bradRecSys = cell(length(bradAreas),1);
    [bradRecSys{1:length(bradAreas)}] = deal('brad');
else
    bradSDF = {}; bradTimes = {}; bradAreas = {}; bradRows = []; bradMonks = {}; bradOtherTasks = {}; bradRecSys = {}; bradSess = {}; bradChans = {};
end

if ismember('tdt',recCrit)
    % Get TDT units
    [outSDF,outSDFTimes,outAreas,outRows,outMonks,outOtherTasks,outSess,outChans,outParams] = klPullKilosortSDFs(task,'-a',areaCrit,'-r',doReward,'-c',clip);
    tdtRecSys = cell(length(outAreas),1);
    [tdtRecSys{1:length(outAreas)}] = deal('tdt');
else
    outSDF = {}; outSDFTimes = {}; outAreas = {}; outRows = []; outMonks = {}; outOtherTasks = {}; tdtRecSys = {}; outSess = {}; outChans = {};
end

if ismember('rich',recCrit)
    % Get Rich's units
    [richSDF,richTimes,richAreas,richRows,richMonks,richOtherTasks,richSess,richChans,richParams] = klPullRichSDFs(task,'-a',areaCrit,'-r',doReward,'-c',clip);
    richRecSys = cell(length(richAreas),1);
    [richRecSys{1:length(richAreas)}] = deal('rich');
else
    richSDF = {}; richTimes = {}; richAreas = {}; richRows = []; richMonks = {}; richOtherTasks = {}; richRecSys = {}; richSess = {}; richChans = {};
end

if ismember('plex',recCrit)
    % Get Plexon units
    [allSDF,allSDFTimes,allAreas,allRows,allMonks,allOtherTasks,allSess,allChans,allParams] = klPullPlexSDFsv2(task,'-a',areaCrit,'-r',doReward,'-c',clip);
    plexRecSys = cell(length(allAreas),1);
    [plexRecSys{1:length(allAreas)}] = deal('plex');
else
    allSDF = {}; allSDFTimes = {}; allAreas = {}; allRows = []; allMonks = {}; allOtherTasks = {}; plexRecSys = {}; allSess = {}; allChans = {};
end



allRecSys = cat(1,plexRecSys,tdtRecSys,richRecSys,bradRecSys);
comboMonks = cat(1,allMonks,outMonks,richMonks,bradMonks);
comboRows = cat(1,allRows,outRows,richRows,bradRows);
comboOtherTasks = cat(1,allOtherTasks,outOtherTasks,richOtherTasks,bradOtherTasks);
% Concatenate
comboSDF = cat(1,allSDF,outSDF,richSDF,bradSDF);
comboSDFTimes = cat(1,allSDFTimes,outSDFTimes,richTimes,bradTimes);
comboAreas = cat(1,allAreas,outAreas,richAreas,bradAreas);
comboSess = cat(1,allSess,outSess,richSess,bradSess);
comboChans = cat(1,allChans,outChans,richChans,bradChans);
comboParams = [allParams,outParams,richParams,bradParams];

switch task
    case {'MG','mg'}
        % Cut out stuff if needed
        badLength = any(cellfun(@(x) ismember(length(x),[0,1]),comboSDFTimes),2);
        comboSDFTimes = comboSDFTimes(~badLength,:);
        comboSDF = comboSDF(~badLength,:);
        comboAreas = comboAreas(~badLength,:);
        comboMonks = comboMonks(~badLength,:);
        allRecSys = allRecSys(~badLength,:);
        comboRows = comboRows(~badLength,:);
        comboOtherTasks = comboOtherTasks(~badLength,:);
        comboSess = comboSess(~badLength,:);
        comboChans = comboChans(~badLength,:);
        comboParams = comboParams(~badLength);
        
        % Align and make a matrix
        allZeros = cellfun(@(x) find(abs(x) == min(abs(x)),1),comboSDFTimes);
        sdfMat = cell(1,size(comboSDF,2)); goodSDFTimes = cell(1,size(comboSDF,2));
        for ii = 1:size(comboSDF,2)
            sdfMat{ii} = klAlignv2(comboSDF(:,ii),allZeros(:,ii));
            goodSDFTimes{ii} = nanmean(klAlignv2(comboSDFTimes(:,ii),allZeros(:,ii)),1);
        end

        allV = sdfMat{1}(:,ismember(goodSDFTimes{1},-200:300));
        allM = sdfMat{2}(:,ismember(goodSDFTimes{2},-300:200));
        catResp = [allV,allM];
        goodRows = sum(isfinite(catResp),2)==size(catResp,2);

        goodSDF = cell(1,length(sdfMat)); 
        for ii = 1:length(sdfMat)
            goodSDF{ii} = sdfMat{ii}(goodRows,:);
        end
        goodAreas = comboAreas(goodRows);
        goodRecSys = allRecSys(goodRows);
        goodRow = comboRows(goodRows);
        goodMonks = comboMonks(goodRows);
        goodOtherTasks = comboOtherTasks(goodRows);
        goodSess = comboSess(goodRows);
        goodChans = comboChans(goodRows);
        goodParams = comboParams(goodRows);
    case {'Cap','Search','Capture','cap','search','capture'}
        fullNans=cellfun(@(x) sum(isnan(x(:)))==numel(x),comboSDF);
        [nanRow,nanCol] = find(fullNans);
        uRows = unique(nanRow);
        goodRows = ones(size(comboSDF,1),1);
        goodRows(uRows) = 0;
        
        goodSDF = comboSDF(logical(goodRows),:);
        goodSDFTimes = comboSDFTimes(logical(goodRows),:);
        goodAreas = comboAreas(logical(goodRows));
        goodRecSys = allRecSys(logical(goodRows));
        goodRow = comboRows(logical(goodRows));
        goodMonks = comboMonks(logical(goodRows));
        goodOtherTasks = comboOtherTasks(logical(goodRows));
        goodSess = comboSess(logical(goodRows));
        goodChans = comboChans(logical(goodRows));
end

