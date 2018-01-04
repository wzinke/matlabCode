function klPreMotorSearch()

close all;

rePull = 0;
whichDistract = 1; %1 = distinguish salient/non-salient distractors, 2 = combine them
doSalient = 0;

if rePull
    [goodSDF,goodSDFTimes,goodAreas] = klPullAllSDFs('cap');
    save('./matlabCode/preMotorAnalysis/fefF2_Search.mat','goodSDF','goodAreas','goodSDFTimes','-v7.3');
else
    load('fefF2_Search.mat');
end

% Normalizing is more difficult here... perhaps a front-end normalization
% is in order but for now (9/11/17) we'll pull the values from the average
% of the conditions

% Get average
avSDF = cellfun(@(x) nanmean(x,1),goodSDF,'UniformOutput',0);
% Make catSDF
[vtCell{1:size(avSDF,1)}] = deal(goodSDFTimes{1});
[mtCell{1:size(avSDF,1)}] = deal(goodSDFTimes{2});
catSDF = cellfun(@(x,y,z,w) [x(:,ismember(y,-200:300)),z(:,ismember(w,-300:200))],...
    avSDF(:,1),vtCell',avSDF(:,2),mtCell','UniformOutput',0);
% Get means and SD
catMeans = cellfun(@nanmean,catSDF,'UniformOutput',0);
catStd = cellfun(@nanstd,catSDF,'UniformOutput',0);
% Make normSDF
normVis = cellfun(@(x,y,z) (x-y)./z,goodSDF(:,1),catMeans,catStd,'UniformOutput',0);
normMov = cellfun(@(x,y,z) (x-y)./z,goodSDF(:,2),catMeans,catStd,'UniformOutput',0);
normSDF = cell(size(goodSDF,1),2);
normSDF(:,1) = normVis;
normSDF(:,2) = normMov;
visTimes = nanmean(cell2mat(goodSDFTimes(:,1)));
movTimes = nanmean(cell2mat(goodSDFTimes(:,2)));

% [normSDF{1},normSDF{2}] = klNormResp(goodSDF{1},goodSDFTimes{1},goodSDF{2},goodSDFTimes{2},'ztr');

% Get a vector of FEF or premotor indices
isFEF = ismember(goodAreas,'FEF');
isF2 = ismember(goodAreas,'F2');

% Get TST for each row (if it exists)
tst = nan(size(normSDF,1));
dst = nan(size(normSDF,1));
for ir = 1:size(normSDF,1)
    switch whichDistract
        case 1
            tst(ir) = getSST(normSDF{ir}(1,:),normSDF{ir}(3,:),visTimes);
            dst(ir) = getSST(normSDF{ir}(3,:),normSDF{ir}(2,:),visTimes);
        case 2
            tst(ir) = getSST(normSDF{ir}(1,:),normSDF{ir}(4,:),visTimes);
    end
end


% Plot TST CDFs for FEF (red) and F2 (blue)
[sortTstFef,~] = sort(tst(isFEF));
[sortTstF2,~] = sort(tst(isF2));
figure(); hold on;
plot(sortTstFef,(1:length(sortTstFef))./sum(isfinite(sortTstFef)),'color',[.8 .2 .2]);
plot(sortTstF2,(1:length(sortTstF2))./sum(isfinite(sortTstF2)),'color',[.2 .2 .8]);

% For those units that have TST, plot T in and D in aligned on stimulus
% FEF on left, F2 on right
% visAlign
switch whichDistract
    case 1
        if doSalient
            fefTargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isFEF & isfinite(tst) & isfinite(dst),1),'UniformOutput',0));
            fefDistIns = cell2mat(cellfun(@(x) x(3,:),normSDF(isFEF & isfinite(tst) & isfinite(dst),1),'UniformOutput',0));
            fefSingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isFEF & isfinite(tst) & isfinite(dst),1),'UniformOutput',0));
            f2TargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isF2 & isfinite(tst) & isfinite(dst),1),'UniformOutput',0));
            f2DistIns = cell2mat(cellfun(@(x) x(3,:),normSDF(isF2 & isfinite(tst) & isfinite(dst),1),'UniformOutput',0));
            f2SingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isF2 & isfinite(tst) & isfinite(dst),1),'UniformOutput',0));
        else
            fefTargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isFEF & isfinite(tst),1),'UniformOutput',0));
            fefDistIns = cell2mat(cellfun(@(x) x(3,:),normSDF(isFEF & isfinite(tst),1),'UniformOutput',0));
            fefSingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isFEF & isfinite(tst),1),'UniformOutput',0));
            f2TargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isF2 & isfinite(tst),1),'UniformOutput',0));
            f2DistIns = cell2mat(cellfun(@(x) x(3,:),normSDF(isF2 & isfinite(tst),1),'UniformOutput',0));
            f2SingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isF2 & isfinite(tst),1),'UniformOutput',0));
        end
    case 2
        fefTargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isFEF & isfinite(tst),1),'UniformOutput',0));
        fefDistIns = cell2mat(cellfun(@(x) x(4,:),normSDF(isFEF & isfinite(tst),1),'UniformOutput',0));
        fefSingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isFEF & isfinite(tst),1),'UniformOutput',0));
        f2TargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isF2 & isfinite(tst),1),'UniformOutput',0));
        f2DistIns = cell2mat(cellfun(@(x) x(4,:),normSDF(isF2 & isfinite(tst),1),'UniformOutput',0));
        f2SingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isF2 & isfinite(tst),1),'UniformOutput',0));
end

figure();
sp(1) = subplot(1,2,1);
pltMeanStd(visTimes,nanmean(fefTargIns,1),nanstd(fefTargIns,[],1)./sqrt(size(fefTargIns,1)),'color',[.8 .2 .2]);
pltMeanStd(visTimes,nanmean(fefDistIns,1),nanstd(fefDistIns,[],1)./sqrt(size(fefDistIns,1)),'color',[.2 .2 .8]);
if whichDistract == 1 && doSalient
    pltMeanStd(visTimes,nanmean(fefSingIns,1),nanstd(fefSingIns,[],1)./sqrt(size(fefSingIns,1)),'color',[.2 .8 .2]);
end
sp(2) = subplot(1,2,2);
pltMeanStd(visTimes,nanmean(f2TargIns,1),nanstd(f2TargIns,[],1)./sqrt(size(f2TargIns,1)),'color',[.8 .2 .2]);
pltMeanStd(visTimes,nanmean(f2DistIns,1),nanstd(f2DistIns,[],1)./sqrt(size(f2DistIns,1)),'color',[.2 .2 .8]);
if whichDistract == 1 && doSalient
    pltMeanStd(visTimes,nanmean(f2SingIns,1),nanstd(f2SingIns,[],1)./sqrt(size(f2SingIns,1)),'color',[.2 .8 .2]);
end
linkaxes(sp,'y');
set(sp,'XLim',[-200,300]);

% Let's plot them individually...
figure();
sp(1) = subplot(1,2,1); hold on;
% plot(visTimes,fefTargIns,'color',[.8 .2 .2]);
plot(visTimes,fefDistIns,'color',[.2 .2 .8]);
plot(visTimes,fefTargIns,'color',[.8 .2 .2]);
if whichDistract == 1 && doSalient
    plot(visTimes,fefSingIns,'color',[.2 .8 .2]);
end
sp(2) = subplot(1,2,2); hold on;
% plot(visTimes,f2TargIns,'color',[.8 .2 .2]);
plot(visTimes,f2DistIns,'color',[.2 .2 .8]);
plot(visTimes,f2TargIns,'color',[.8 .2 .2]);
if whichDistract == 1 && doSalient
    plot(visTimes,f2SingIns,'color',[.2 .8 .2]);
end
linkaxes(sp,'y');
set(sp,'XLim',[-200,300]);

% Let's try another where we plot the difference functions
figure();
sp(1) = subplot(2,2,1); hold on;
pltMeanStd(visTimes,nanmean(fefTargIns-fefDistIns,1),nanstd(fefTargIns-fefDistIns,[],1)./sqrt(size(fefTargIns,1)),'color',[.8 .2 .2]);
pltMeanStd(visTimes,nanmean(f2TargIns-f2DistIns,1),nanstd(f2TargIns-f2DistIns,[],1)./sqrt(size(f2TargIns,1)),'color',[.2 .2 .8]);
sp(2) = subplot(2,2,2); hold on;
plot(visTimes,fefTargIns-fefDistIns,'color',[.8 .2 .2]);
plot(visTimes,f2TargIns-f2DistIns,'color',[.2 .2 .8]);
set(sp,'XLim',[-200 300]);



% For those units that have TST, plot T in and D in aligned on saccade
% FEF on left, F2 on right
switch whichDistract
    case 1
        if doSalient
            fefTargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isFEF & isfinite(tst) & isfinite(dst),2),'UniformOutput',0));
            fefDistIns = cell2mat(cellfun(@(x) x(3,:),normSDF(isFEF & isfinite(tst) & isfinite(dst),2),'UniformOutput',0));
            fefSingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isFEF & isfinite(tst) & isfinite(dst),2),'UniformOutput',0));
            f2TargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isF2 & isfinite(tst) & isfinite(dst),2),'UniformOutput',0));
            f2DistIns = cell2mat(cellfun(@(x) x(3,:),normSDF(isF2 & isfinite(tst) & isfinite(dst),2),'UniformOutput',0));
            f2SingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isF2 & isfinite(tst) & isfinite(dst),2),'UniformOutput',0));
        else
            fefTargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isFEF & isfinite(tst),2),'UniformOutput',0));
            fefDistIns = cell2mat(cellfun(@(x) x(3,:),normSDF(isFEF & isfinite(tst),2),'UniformOutput',0));
            fefSingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isFEF & isfinite(tst),2),'UniformOutput',0));
            f2TargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isF2 & isfinite(tst),2),'UniformOutput',0));
            f2DistIns = cell2mat(cellfun(@(x) x(3,:),normSDF(isF2 & isfinite(tst),2),'UniformOutput',0));
            f2SingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isF2 & isfinite(tst),2),'UniformOutput',0));
        end
    case 2
        fefTargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isFEF & isfinite(tst),2),'UniformOutput',0));
        fefDistIns = cell2mat(cellfun(@(x) x(4,:),normSDF(isFEF & isfinite(tst),2),'UniformOutput',0));
        fefSingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isFEF & isfinite(tst),2),'UniformOutput',0));
        f2TargIns = cell2mat(cellfun(@(x) x(1,:),normSDF(isF2 & isfinite(tst),2),'UniformOutput',0));
        f2DistIns = cell2mat(cellfun(@(x) x(4,:),normSDF(isF2 & isfinite(tst),2),'UniformOutput',0));
        f2SingIns = cell2mat(cellfun(@(x) x(2,:),normSDF(isF2 & isfinite(tst),2),'UniformOutput',0));
end

figure();
sp(1) = subplot(1,2,1);
pltMeanStd(movTimes,nanmean(fefTargIns,1),nanstd(fefTargIns,[],1)./sqrt(size(fefTargIns,1)),'color',[.8 .2 .2]);
pltMeanStd(movTimes,nanmean(fefDistIns,1),nanstd(fefDistIns,[],1)./sqrt(size(fefDistIns,1)),'color',[.2 .2 .8]);
if whichDistract == 1 && doSalient
    pltMeanStd(movTimes,nanmean(fefSingIns,1),nanstd(fefSingIns,[],1)./sqrt(size(fefSingIns,1)),'color',[.2 .8 .2]);
end
sp(2) = subplot(1,2,2);
pltMeanStd(movTimes,nanmean(f2TargIns,1),nanstd(f2TargIns,[],1)./sqrt(size(f2TargIns,1)),'color',[.8 .2 .2]);
pltMeanStd(movTimes,nanmean(f2DistIns,1),nanstd(f2DistIns,[],1)./sqrt(size(f2DistIns,1)),'color',[.2 .2 .8]);
if whichDistract == 1 && doSalient
    pltMeanStd(movTimes,nanmean(f2SingIns,1),nanstd(f2SingIns,[],1)./sqrt(size(f2SingIns,1)),'color',[.2 .8 .2]);
end
linkaxes(sp,'y');
set(sp,'XLim',[-300,300]);

% Let's plot them individually...
figure();
sp(1) = subplot(1,2,1); hold on;
plot(movTimes,fefTargIns,'color',[.8 .2 .2]);
plot(movTimes,fefDistIns,'color',[.2 .2 .8]);
if whichDistract == 1 && doSalient
    plot(movTimes,fefSingIns,'color',[.2 .8 .2]);
end
sp(2) = subplot(1,2,2); hold on;
plot(movTimes,f2TargIns,'color',[.8 .2 .2]);
plot(movTimes,f2DistIns,'color',[.2 .2 .8]);
if whichDistract == 1 && doSalient
    plot(movTimes,f2SingIns,'color',[.2 .8 .2]);
end
linkaxes(sp,'y');
set(sp,'XLim',[-300,300]);

% Let's try another where we plot the difference functions
figure();
sp(1) = subplot(2,2,1); hold on;
pltMeanStd(movTimes,nanmean(fefTargIns-fefDistIns,1),nanstd(fefTargIns-fefDistIns,[],1)./sqrt(size(fefTargIns,1)),'color',[.8 .2 .2]);
pltMeanStd(movTimes,nanmean(f2TargIns-f2DistIns,1),nanstd(f2TargIns-f2DistIns,[],1)./sqrt(size(f2TargIns,1)),'color',[.2 .2 .8]);
sp(2) = subplot(2,2,2); hold on;
plot(movTimes,fefTargIns-fefDistIns,'color',[.8 .2 .2]);
plot(movTimes,f2TargIns-f2DistIns,'color',[.2 .2 .8]);
set(sp,'XLim',[-300 300]);


end


function [sst, diffSDF] = getSST(in,across,times)
    % Get difference and z-score
    diffSDF = in - across;
    blMean = nanmean(diffSDF(times >= -50 & times <= 50));
    blStd = nanstd(diffSDF(times >= -50 & times <= 50));
    zSDF = (diffSDF-blMean)./blStd;
    
    % Get the first time that zSDF >= 2 for more than 15ms
    overFive = klGetConsecutive(zSDF >= 5);
    overTwo = klGetConsecutive(zSDF >= 2);
    firstOver = find(overFive >= 15 & times >= 0,1);
    if overTwo(firstOver) >= 15
        sst = times(find(overTwo(1:firstOver) == 0,1,'last'));
    else
        sst = nan;
    end    
end
