%% Plots means with filled SDs
%  In development - will add varargin and options later

function [line, shape] = pltMeanStd(xVals,means,stds,varargin)

% Set defaults
alpha = .5;
lineWidth = 1;
color = 'k';
edge = 0;
doFill = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'color'}
            color = varargin{varStrInd(iv)+1};
        case {'width','linewidth'}
            lineWidth = varargin{varStrInd(iv)+1};
        case {'alpha','facealpha'}
            alpha = varargin{varStrInd(iv)+1};
        case {'edge'},
            edge = varargin{varStrInd(iv)+1};
        case {'-f'}
            doFill = varargin{varStrInd(iv)+1};
    end
end

if edge,
    eColor = color;
else
    eColor = 'none';
end

if size(stds,2) == 2 && size(stds,1) > 2,
    stds = stds';
end

if size(stds,1) == 1,
    upMeans = means+stds;
    downMeans = means-stds;
elseif size(stds,1) == 2,
    upMeans = stds(1,:);%+means
    downMeans = stds(2,:);%+means
end

if size(stds,1) == 1,
    if sum(isnan(means)) == length(means) || sum(isnan(stds)) == length(stds),
        means(length(means):max([length(means),length(stds)])) = nan;
        stds(length(stds):max([length(means),length(stds)])) = nan;
    end
end

if any(isnan(upMeans)) || any(isnan(downMeans))
    warning('NaN detected in SD values...\n\tFill may not work properly\n');
    upMeans(isnan(stds) | isnan(means))    = [];
    downMeans(isnan(stds) | isnan(means))  = [];
    xVals(isnan(stds) | isnan(means))      = [];
    means(isnan(stds) | isnan(means))      = [];
    %keyboard;
end

xValsFill = [xVals,xVals(end:-1:1)];
stdFill = [upMeans,downMeans(end:-1:1)];

hold on;
if doFill
    shape = fill(xValsFill,stdFill,color,'FaceAlpha',alpha,'edgecolor',eColor','edgealpha',alpha); hold on;
else
    plot(xVals,upMeans,'color',color,'linewidth',lineWidth);
    plot(xVals,downMeans,'color',color,'linewidth',lineWidth);
    plot([xVals(1),xVals(1)],[upMeans(1),downMeans(1)],'color',color,'linewidth',lineWidth);
    plot([xVals(end),xVals(end)],[upMeans(end),downMeans(end)],'color',color,'linewidth',lineWidth);
end 
    
line = plot(xVals,means,'color',color,'linewidth',lineWidth);
