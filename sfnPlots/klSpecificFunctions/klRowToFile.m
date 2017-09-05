function [outPath, outFile] = klRowToFile(row,varargin)

% Set defaults
file = 'klDataBookKeeping_mg.xlsx';
monk = 'Gauss';
task = 'MG';
path = 'Y:/Users/Wolf/ephys_db';
fType = 'DSP';
wave = 0;
doLFP = 0;
reload = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'file','-f'}
            file = varargin{varStrInd(iv)+1};
        case {'monk','-m'}
            monk = varargin{varStrInd(iv)+1};
        case {'task','-t'}
            task = varargin{varStrInd(iv)+1};
        case {'ftype'}
            fType = varargin{varStrInd(iv)+1};
        case {'wave','-w'}
            wave  = varargin{varStrInd(iv)+1};
            if wave, doLFP = 0; end
        case {'-r','reload'}
            reload = varargin{varStrInd(iv)+1};
        case {'-l','lfp'},
            doLFP = varargin{varStrInd(iv)+1};
            if doLFP, wave = 0; end
    end
end

global excelNum excelAll

if isempty(excelNum) || isempty(excelAll) || reload
    [excelNum,~,excelAll] = xlsread(file,monk);
end
hRow = find(strcmp(excelAll(:,2),'Name'),1);
sessCol = find(strcmp(excelAll(hRow,:),'Name'),1);
chanCol = find(strcmp(excelAll(hRow,:),'chanCode'),1);

chanCode = excelAll{row,chanCol};
if strcmp(fType,'LFP'),
    chanCode = ['LFP',chanCode(4:5)];
end
outPath = sprintf('%s/%s/%s/%s/%s/',path,monk,excelAll{row,sessCol},fType,chanCode);

if wave,
    outPath = [outPath,'waves/'];
end

outFileStruct = dir(sprintf('%s/%s_%s_*%s*',outPath,excelAll{row,sessCol},chanCode,task));
for ifile = 1:length(outFileStruct)
    outFile{ifile} = outFileStruct(ifile).name(1:(end-4));
end

nTries = 0;
% Try up to 20 times to get the file
while ~exist('outFile','var') && nTries < 20,
    outPath = sprintf('%s/%s/%s/%s/%s/',path,monk,excelAll{row,sessCol},fType,excelAll{row,chanCol});
    if wave,
        outPath = [outPath,'waves/'];
    end

    outFileStruct = dir(sprintf('%s/%s_%s_*%s*',outPath,excelAll{row,sessCol},chanCode,task));
    for ifile = 1:length(outFileStruct)
        outFile{ifile} = outFileStruct(ifile).name(1:(end-4));
    end
    nTries = nTries + 1;
end

if ~exist('outFile','var'),
%     keyboard
    outFile{1} = nan;
end