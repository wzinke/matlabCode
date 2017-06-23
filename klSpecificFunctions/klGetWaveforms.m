function grpWaves = klGetWaveforms(xlRows,varargin)

% Set defaults
monk = 'Gauss';
task = 'MG';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'monk','-m'}
            monk = varargin{varStrInd(iv)+1};
    end
end

if sum(ismember(xlRows,[0,1])) == length(xlRows),
    xlRows = find(xlRows);
end

% Set constants
xlFile = 'klDataBookKeeping_mg.xlsx';

% Load in excel file
global excelNum excelAll
if isempty(excelNum) || isempty(excelAll)
    [excelNum,~,excelAll] = xlsread(xlFile,monk);
end

% Start xl row loop
grpWaves = nan(length(xlRows),32);
for ir = 1:length(xlRows),
    [path,file] = klRowToFile(xlRows(ir),'-w',1,'-m',monk,'-t',task);
    if isempty(file), continue; end;
    load([path,file{1}]);
    grpWaves(ir,:) = nanmean(wave.waves,1);
end