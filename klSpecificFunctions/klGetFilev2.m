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
% load('sortedWaves.mat');
load allSorts.mat

% Start xl row loop
grpWaves = cell(length(xlRows),1);
grpWvTimes = cell(length(xlRows),1);
grpSDF = cell(length(xlRows),2);
grpStd = cell(length(xlRows),2);
grpSDFTimes = cell(length(xlRows),2);
normMat = nan(length(xlRows),2);
for ir = 1:length(xlRows),
    [path,file] = klRowToFile(xlRows(ir),'-m',monk,'-t',task);
    if isempty(file), continue; end;
    load([path,file{1}]);
    saccTime = Task.GoCue+Task.SRT;
    [tmpSDF,tmpTimes] = klSpkRatev2(spiketimes,'-q',1);
    bl = nanmean(nanmean(tmpSDF(:,ismember(tmpTimes,blWind)),zDim(1)),zDim(2));
    [rf antirf]= getRFv2(spiketimes,Task,'-b',bl);
    vRF = rf(1); mRF = rf(2);
    vARF = antirf(1); mARF = antirf(2);
    if doRF,
        vCrit = Task.TargetLoc == vRF;
        mCrit = Task.TargetLoc == mRF;
        vaCrit = Task.TargetLoc == vARF;
        maCrit = Task.TargetLoc == mARF;
    else
        vCrit = true(size(spiketimes,1),1);
        mCrit = true(size(spiketimes,1),1);
    end
    
    [rawVisSDF,grpSDFTimes{ir,3}] = klSpkRatev2(spiketimes(vCrit,:),'-q',1);
    [mSDF, grpSDFTimes{ir,2}] = klSpkRatev2(spiketimes(mCrit,:)-repmat(Task.SRT(mCrit,:)+Task.GoCue(mCrit,:),1,size(spiketimes,2)),'-q',1);
    if clip,
        spiketimes(spiketimes > repmat(saccTime,1,size(spiketimes,2))) = nan;
    end
    [vSDF, grpSDFTimes{ir,1}] = klSpkRatev2(spiketimes(vCrit,:),'-q',1);
    normMat(ir,1) = nanmean(nanmean(vSDF(:,ismember(grpSDFTimes{ir,1},blWind)),zDim(1)),zDim(2));
    switch zType,
        case 'baseline',
            normMat(ir,2) = nanstd(nanmean(rawVisSDF(:,ismember(grpSDFTimes{ir,1},blWind)),zDim(1)),[],zDim(2));
        case 'trial',
            normMat(ir,2) = nanstd(nanmean(rawVisSDF,zDim(1)),[],zDim(2));
    end
    normMat(ir,3) = nanmax(abs(nanmean(vSDF(:,ismember(grpSDFTimes{ir,1},vWind)),1)-normMat(ir,1)));
    normMat(ir,4) = nanmax(abs(nanmean(mSDF(:,ismember(grpSDFTimes{ir,2},mWind)),1)-normMat(ir,1)));
    switch groupVar
        case 'none'
            vTrGroup = ones(size(vSDF,1),1);
            mTrGroup = ones(size(mSDF,1),1);
        case {'loc'}
            vTrGroup = Task.TargetLoc;
            mTrGroup = Task.TargetLoc;
        case {'targ','target','tst'}
            vTrGroup = Task.TargetLoc == vRF;
            mTrGroup = Task.TargetLoc == mRF;
    end
    
    uGroup = unique(vTrGroup);
    for ig = 1:length(uGroup)
        grpSDF{ir,1}(ig,:) = nanmean(vSDF(vTrGroup == uGroup(ig),:),1);
        grpSDF{ir,3}(ig,:) = nanmean(rawVisSDF(vTrGroup == uGroup(ig),:),1);
        grpStd{ir,1}(ig,:) = nanstd(vSDF(vTrGroup == uGroup(ig),:),[],1)./sqrt(sum(isfinite(vSDF(vTrGroup == uGroup(ig),:)),1));
        grpStd{ir,3}(ig,:) = nanstd(rawVisSDF(vTrGroup == uGroup(ig),:),[],1)./sqrt(sum(isfinite(rawVisSDF(vTrGroup == uGroup(ig),:)),1));
    end
    uGroup = unique(mTrGroup);
    for ig = 1:length(uGroup)
        grpSDF{ir,2}(ig,:) = nanmean(mSDF(mTrGroup == uGroup(ig),:),1);
        grpStd{ir,2}(ig,:) = nanstd(mSDF(mTrGroup == uGroup(ig),:),[],1)./sqrt(sum(isfinite(mSDF(mTrGroup == uGroup(ig),:)),1));
    end
    
   [path,file] = klRowToFile(xlRows(ir),'-w',1,'-m',monk,'-t',task);
   if isempty(file), continue; end;
   load([path,file{1}]);
%    grpWaves(ir,:) = nanmean(wave.waves,1);
grpWaves{ir} = sortWaves.(monk).(sprintf('s%s',excelAll{xlRows(ir),2}(~ismember(excelAll{xlRows(ir),2},'-')))).(excelAll{xlRows(ir),5}).apWave;
grpWvTimes{ir} = sortWaves.(monk).(sprintf('s%s',excelAll{xlRows(ir),2}(~ismember(excelAll{xlRows(ir),2},'-')))).(excelAll{xlRows(ir),5}).apTimes;
end