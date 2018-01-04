clear all; close all;

rePull = 0;

if rePull
    [goodSDF,goodSDFTimes,goodAreas] = klPullAllSDFs('mg');
else
    load('fefF2MG.mat');
end

% Classify them
[numTradTypes] = klTradClassify(goodSDF,goodSDFTimes);

% Use ZTR normalization
[normSDF{1},normSDF{2}] = klNormResp(goodSDF{1},goodSDFTimes{1},goodSDF{2},goodSDFTimes{2},'ztr');

% Figure 1-4: Plot Vis, Vismov, Mov, None Cells: 
% FEF on top, F2 below, Vis left, Mov right
for ii = 1:4
    figure(ii);
    sp(ii,1) = subplot(2,2,1);
    plot(goodSDFTimes{1},normSDF{1}(numTradTypes == ii & ismember(goodAreas,'FEF'),:),'color',[.3 .3 .3]);
    set(gca,'XLim',[-200 300]);
    sp(ii,2) = subplot(2,2,2);
    plot(goodSDFTimes{2},normSDF{2}(numTradTypes == ii  & ismember(goodAreas,'FEF'),:),'color',[.3 .3 .3]);
    set(gca,'XLim',[-300 200]);
    sp(ii,3) = subplot(2,2,3);
    plot(goodSDFTimes{1},normSDF{1}(numTradTypes == ii & ismember(goodAreas,'F2'),:),'color',[.3 .3 .3]);
    set(gca,'XLim',[-200 300]);
    sp(ii,4) = subplot(2,2,4);
    plot(goodSDFTimes{2},normSDF{2}(numTradTypes == ii  & ismember(goodAreas,'F2'),:),'color',[.3 .3 .3]);
    set(gca,'XLim',[-300 200]);
    
end
linkaxes(sp,'y');
clear sp;


% Figure 100 + 1-4: Plot Mean Vis, Vismov, Mov, None Cells: 
% FEF on top, F2 below, Vis left, Mov right
for ii = 1:4
    figure(ii+100);
    sp(ii,1) = subplot(2,2,1);
    plot(goodSDFTimes{1},nanmean(normSDF{1}(numTradTypes == ii & ismember(goodAreas,'FEF'),:),1),'color',[.3 .3 .3]);
    set(gca,'XLim',[-200 300]);
    sp(ii,2) = subplot(2,2,2);
    plot(goodSDFTimes{2},nanmean(normSDF{2}(numTradTypes == ii  & ismember(goodAreas,'FEF'),:),1),'color',[.3 .3 .3]);
    set(gca,'XLim',[-300 200]);
    sp(ii,3) = subplot(2,2,3);
    plot(goodSDFTimes{1},nanmean(normSDF{1}(numTradTypes == ii & ismember(goodAreas,'F2'),:),1),'color',[.3 .3 .3]);
    set(gca,'XLim',[-200 300]);
    sp(ii,4) = subplot(2,2,4);
    plot(goodSDFTimes{2},nanmean(normSDF{2}(numTradTypes == ii  & ismember(goodAreas,'F2'),:),1),'color',[.3 .3 .3]);
    set(gca,'XLim',[-300 200]);
    
end

linkaxes(sp,'y');
clear sp;