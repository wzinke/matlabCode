function [combGrp, clipTimes] = klPlotGrp(rows,varargin)

% Set defaults
monk    = 'Gauss';
alignOn = {'StimOnset','SRT'};
split   = [];
blWind  = [-250,0];
tWind   = -200:1300;
plotChans = 1;
blDim = 2;
clip = 1;
normType = 'z';
%figStart = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'monk','-m'}
            monk = varargin{varStrInd(iv)+1};
        case {'plot','-p'}
            plotChans = varargin{varStrInd(iv)+1};
        case {'fig','-f'}
            figStart = varargin{varStrInd(iv)+1};
        case {'align','-a'}
            alignOn = varargin{varStrInd(iv)+1};
        case {'bDim','-b'},
            blDim = varargin{varStrInd(iv)+1};
        case {'clip','-c'},
            clip = varargin{varStrInd(iv)+1};
        case {'-n','norm'},
            normType = varargin{varStrInd(iv)+1};
    end
end

% Check input
if isempty(rows) || sum(rows) == 0,
    combGrp = {};
    clipTimes = {};
    return;
end

% Set constants
xlFile = 'klDataBookKeeping_mg.xlsx';
dims  = [1,2];

% Read in excel file
global excelNum excelAll     % Make global for downstream functions
[excelNum,~,excelAll] = xlsread(xlFile,monk); %Reread for this function because monk might have changed

% Loop through rows
for ir = 1:length(rows)
    % Load this row data
    [path,file] = klRowToFile(rows(ir),'-m',monk);
    fStruct = load([path,file{1}]);
    Task = fStruct.Task;
    spiketimes = fStruct.spiketimes;
    
    % Get baseline mean/std
    [sdfRaw,sdfRawT]    = klSpkRatev2(spiketimes,'-q',1);
    blInds = sdfRawT >= blWind(1) & sdfRawT <= blWind(2);
    blSpks = nanmean(sdfRaw(:,blInds),blDim);
    blMean = nanmean(blSpks); %blMean = nanmean(blMean,dims(~ismember(dims,blDim)));
    blStd  = nanstd(blSpks);
    
    % Get mean SDF
    for ia = 1:length(alignOn)
        alignTimes          = Task.(alignOn{ia}); if strcmp(alignOn{ia},'SRT'),alignTimes = alignTimes + Task.GoCue; end
        if sum(isnan(alignTimes)) == numel(alignTimes),
            if strcmp(alignOn{ia},'SRT'),
                alignTimes = Task.SaccEnd;
                if sum(isnan(alignTimes)) == numel(alignTimes),
                    keyboard
                end
            else
                keyboard
            end
        end
        
        clipSpks = spiketimes;
        if clip,
            if ia == 1,
                clipEv = Task.(alignOn{2});
                clipSpks(spiketimes > repmat(clipEv,1,size(spiketimes,2))) = nan;
            elseif ia == 2,
                clipEv = Task.(alignOn{1});
                clipSpks(spiketimes < repmat(clipEv,1,size(spiketimes,2))) = nan;
            end
        end
        
        [sdf,sdfTimes]      = klSpkRatev2(clipSpks-repmat(alignTimes,1,size(spiketimes,2)),'-q',1);
        if isempty(split)
            switch normType
                case {'-z','z'}
                    normSDF = normData(sdf,'-z',[blMean,blStd]);
                case {'max'}
                    mnSDF = nanmean(sdf,1);
                    maxVal = max(abs(mnSDF(:,ismember(sdfTimes,tWind))-blMean));
                    normSDF = (mnSDF-blMean)./maxVal;
            end
            grpSDF{ir,ia}       = nanmean(normSDF,1);
            grpSDFTimes{ir,ia}  = sdfTimes;
        end
        
    end
    
end

% Get the SDF lengths for the included cells
grpLens = cellfun(@length,grpSDF);
minLens = min(grpLens,[],1);
for ia = 1:length(alignOn),
    [minCell{1:size(grpLens,1),1}] = deal(minLens(ia));
    [dimCell{1:size(grpLens,1),1}] = deal(ia);
    clipGrp(:,ia) = cellfun(@clipToMin,grpSDF(:,ia),minCell,dimCell,'UniformOutput',0);
    clipTimes(:,ia) = cellfun(@clipToMin,grpSDFTimes(:,ia),minCell,dimCell,'UniformOutput',0);
    combGrp{ia} = cell2mat(clipGrp(:,ia));
    
    if plotChans
        if exist('figStart','var'),
            figure(figStart+(ia-1));
        end
        if size(combGrp{ia},1) < 10,
            hold on;
            plot(nanmean(cell2mat(clipTimes(:,ia)),1),combGrp{ia});
        end
        pltMeanStd(nanmean(cell2mat(clipTimes(:,ia)),1),nanmean(combGrp{ia},1),nanstd(combGrp{ia},1)./sqrt(size(combGrp{ia},1)),'k');
        hline(2);
    end
    
end
% keyboard
function clipped = clipToMin(raw,minVal,dim)
    % Dim 1 is forward, dim 2 = backwards
    if dim == 1,
        clipped = raw(1:minVal);
    elseif dim == 2,
        clipped = raw((end-(minVal-1)):end);
    end
end

end