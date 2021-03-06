function [grpSDF, grpSDFTimes, grpWaves, grpWvTimes, normMat, grpStd] = klGetFilev2(xlRows,varargin)

% Set defaults
monk = 'Gauss';
groupVar = 'none';
task = 'MG';
clip = 1;
blWind = -200:-100;
vWind = -100:400;
mWind = -400:100;
zDim = [1,2];
zType = 'baseline';
doRF = 0;
doDiff = 0;
diffType = 'opposite';
plexWvTimes = ((1:32)-9).*25;
tdtWvTimes = ((1:32)-9).*(1000/24414).*1000;
splineTimes = spline(1:32,plexWvTimes,1:.1:32);
print = 1;
nPrint = 20;
recSys = 'plex';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'monk','-m'}
            monk = varargin{varStrInd(iv)+1};
        case {'task','-t'}
            task = varargin{varStrInd(iv)+1};
        case {'rf','-rf'}
            doRF = varargin{varStrInd(iv)+1};
        case {'-g','group'},
            groupVar = varargin{varStrInd(iv)+1};
        case {'-z'}
            zDim = varargin{varStrInd(iv)+1};
            if length(zDim) == 1,
                zDirs = [1,2];
                zDim(2) = zDim(1);
                zDim(1) = zDirs(~ismember(zDirs,zDim));
            end
        case {'ztype'}
            zType = varargin{varStrInd(iv)+1};
        case {'-p'},
            nPrint = varargin{varStrInd(iv)+1};
            print = 1;
        case {'sys'},
            recSys = varargin{varStrInd(iv)+1};
        case {'diff'},
            doDiff = varargin{varStrInd(iv)+1};
    end
end

if sum(ismember(xlRows,[0,1])) == length(xlRows),
    xlRows = find(xlRows);
end

% Set constants
xlFile = 'klDataBookKeeping_mg.xlsx';
tdtFile = 'klTDTBookKeeping.xlsx';

% Load in excel file
global excelNum excelAll masterSorts
if isempty(excelNum) || isempty(excelAll)
    switch recSys
        case {'plex','p','-p'},
            [excelNum,~,excelAll] = xlsread(xlFile,monk);
        case {'tdt','t','-t'},
            [excelNum,~,excelAll] = xlsread(tdtFile,monk);
    end
end
% load('sortedWaves.mat');
% load allSorts.mat
load meanSorts.mat
load masterSorts

% Start xl row loop
grpWaves = cell(length(xlRows),1);
grpWvTimes = cell(length(xlRows),1);
grpSDF = cell(length(xlRows),2);
grpStd = cell(length(xlRows),2);
grpSDFTimes = cell(length(xlRows),2);
normMat = nan(length(xlRows),2);
for ir = 1:length(xlRows),
    if print && mod(ir,nPrint) == 0,
        fprintf('Working on row %d (of %d)...\n',ir,length(xlRows));
    end
    
    switch recSys
        case {'plex','p','-p'},
            [path,file] = klRowToFile(xlRows(ir),'-m',monk,'-t',task);
            load([path,file{1}]);
        case {'tdt','t','-t'},
            [path,file] = tdtRowToFile(xlRows(ir),'-m',monk,'-r',1);
            load([path,file{1}]);
            [path,file] = tdtRowToFile(xlRows(ir),'-m',monk,'-t',1);
            load([path,file{1}]);
            spiketimes = spikes.spiketimes;
    end
    
    if isempty(file), continue; end;
    Task=klCutTask(Task,ismember(Task.TaskType,'MG'));
    spiketimes = spiketimes(ismember(Task.TaskType,'MG'),:);
    
    if sum(isnan(spiketimes(:))) == numel(spiketimes), 
        continue; 
    end
    
    saccTime = Task.GoCue+Task.SRT;
    [tmpSDF,tmpTimes] = klSpkRatev2(spiketimes,'-q',1);
    bl = nanmean(nanmean(tmpSDF(:,ismember(tmpTimes,blWind)),zDim(1)),zDim(2));
    [rf, antirf]= getRFv2(spiketimes,Task,'-b',bl);
    vRF = rf(1); mRF = rf(2);
    vARF = antirf(1); mARF = antirf(2);
    
    if doRF || doDiff,
        vCrit = Task.TargetLoc == vRF;
        mCrit = Task.TargetLoc == mRF;
        switch diffType,
            case {'arf'},
                vaCrit = Task.TargetLoc == vARF;
                maCrit = Task.TargetLoc == mARF;
            case {'opposite'},
                vaCrit = Task.TargetLoc == mod(vRF+180,360);
                maCrit = Task.TargetLoc == mod(mRF+180,360);
        end
    else
        vCrit = true(size(spiketimes,1),1);
        mCrit = true(size(spiketimes,1),1);
    end
    
    [rawVisSDF,grpSDFTimes{ir,3}] = klSpkRatev2(spiketimes,'-q',1);
    [rawMovSDF,grpSDFTimes{ir,2}] = klSpkRatev2(spiketimes-repmat(Task.SRT+Task.GoCue,1,size(spiketimes,2)),'-q',1);
    if clip,
        spiketimes(spiketimes > repmat(saccTime,1,size(spiketimes,2))) = nan;
    end
    [clipVisSDF,grpSDFTimes{ir,1}] = klSpkRatev2(spiketimes,'-q',1);
    
    normMat(ir,1) = nanmean(nanmean(clipVisSDF(:,ismember(grpSDFTimes{ir,1},blWind)),zDim(1)),zDim(2));
    switch zType,
        case 'baseline',
            normMat(ir,2) = nanstd(nanmean(rawVisSDF(:,ismember(grpSDFTimes{ir,1},blWind)),zDim(1)),[],zDim(2));
        case 'trial',
            if sum(vCrit),
                normMat(ir,2) = nanstd(nanmean(rawVisSDF(vCrit,:),zDim(1)),[],zDim(2));
            else
                normMat(ir,2) = nan;
            end
    end
    if sum(vCrit),
        normMat(ir,3) = nanmax(abs(nanmean(clipVisSDF(:,ismember(grpSDFTimes{ir,1},vWind)),1)-normMat(ir,1)));
        normMat(ir,4) = nanmax(abs(nanmean(rawMovSDF(:,ismember(grpSDFTimes{ir,2},mWind)),1)-normMat(ir,1)));
    else
        normMat(ir,3:4) = [nan,nan];
    end
    
    if doDiff,
        grpSDF{ir,1} = nanmean(clipVisSDF(vCrit,:),1) - nanmean(clipVisSDF(vaCrit,:),1);
        grpSDF{ir,2} = nanmean(rawMovSDF(mCrit,:),1) - nanmean(rawMovSDF(maCrit,:),1);
        grpSDF{ir,3} = nanmean(rawVisSDF(vCrit,:),1) - nanmean(rawVisSDF(vaCrit,:),1);
    else
        grpSDF{ir,1} = nanmean(clipVisSDF(vCrit,:),1);
        grpSDF{ir,2} = nanmean(rawMovSDF(mCrit,:),1);
        grpSDF{ir,3} = nanmean(rawVisSDF(vCrit,:),1);
    end
    
    if isfield(masterSorts,monk) && ~ismember(recSys,{'tdt','t','-t'}),
        grpWaves{ir} = masterSorts.(monk)(xlRows(ir)).wave;
        grpWvTimes{ir} = masterSorts.(monk)(xlRows(ir)).time;
    else
        switch recSys
            case {'plex','p','-p'},
                [path,file] = klRowToFile(xlRows(ir),'-m',monk,'-t',task,'-w',1);
                load([path,file{1}]);
                grpWaves{ir} = spline(plexWvTimes,nanmean(wave.waves,1),splineTimes);
                grpWvTimes{ir} = splineTimes;
            case {'tdt','t','-t'},
                tmpSplWaves = spline(tdtWvTimes,nanmean(spikes.waves,1),splineTimes);
                tmpSplTimes = splineTimes;
                [grpWaves{ir},grpWvTimes{ir}] = klTroughAlignv5(tmpSplWaves,tmpSplTimes,0);
        end
    end
end