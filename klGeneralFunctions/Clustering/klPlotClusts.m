function klPlotClusts(obs,obsX,clusts,tmpK,varargin)

% Set defaults
slow = 0;
indiv = 0;
doSD = 1;
doMedian = 0;
doSeparate = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-s','slow'}
            slow = varargin{varStrInd(iv)+1};
        case {'-i','indiv'},
            indiv = varargin{varStrInd(iv)+1};
        case {'-x','xlim'},
            xlim = varargin{varStrInd(iv)+1};
        case {'-sd'},
            doSD = varargin{varStrInd(iv)+1};
        case {'-c'},
            colors = varargin{varStrInd(iv)+1};
        case {'-m'},
            doMedian = varargin{varStrInd(iv)+1};
        case {'-sep','separate'},
            doSeparate = varargin{varStrInd(iv)+1};
    end
end
            
% Get colors
if ~exist('colors','var') || (size(colors,1) ~= tmpK && ~ischar(colors)) || (length(colors) ~= tmpK && ischar(colors)),
    if tmpK > 5,
        colors = jet(tmpK);
    else
        colors = [
            .8 .2 .2
            .2 .8 .2
            .2 .2 .8
            .2 .8 .8
            .8 .2 .8
            ];
    end
end

% Get Axes handle
axes(gca);

% Plot in the loop
for ik = 1:tmpK,
    if sum(clusts(:,tmpK)==ik),
        if indiv,
            if doSeparate,
                figure(1000+ik);
            end
            plot(obsX,obs(clusts(:,tmpK)==ik,:),'color',colors(ik,:)); hold on;
        elseif ~indiv && ~doSD
            if doMedian,
                plot(obsX,nanmedian(obs(clusts(:,tmpK)==ik,:),1),'color',colors(ik,:)); hold on;
            else
                plot(obsX,nanmean(obs(clusts(:,tmpK)==ik,:),1),'color',colors(ik,:)); hold on;
            end
        else
            if doMedian,
                pltMeanStd(obsX,nanmedian(obs(clusts(:,tmpK)==ik,:),1),nanstd(obs(clusts(:,tmpK)==ik,:),1)./sqrt(sum(clusts(:,tmpK)==ik)),'color',colors(ik,:)); hold on;
            else
                pltMeanStd(obsX,nanmean(obs(clusts(:,tmpK)==ik,:),1),nanstd(obs(clusts(:,tmpK)==ik,:),1)./sqrt(sum(clusts(:,tmpK)==ik)),'color',colors(ik,:)); hold on;
            end
%             errorbar(obsX,nanmean(obs(clusts(:,tmpK)==ik,:),1),nanstd(obs(clusts(:,tmpK)==ik,:),1)./sqrt(sum(clusts(:,tmpK)==ik)),'color',colors(ik,:)); hold on;
%             plot(obsX,nanmean(obs(clusts(:,tmpK)==ik,:),1),'color',colors(ik,:)); hold on;
%             plot(obsX,nanmean(obs(clusts(:,tmpK)==ik,:),1)+nanstd(obs(clusts(:,tmpK)==ik,:),1)./sqrt(sum(clusts(:,tmpK)==ik)),'color',colors(ik,:)); hold on;
%             plot(obsX,nanmean(obs(clusts(:,tmpK)==ik,:),1)-nanstd(obs(clusts(:,tmpK)==ik,:),1)./sqrt(sum(clusts(:,tmpK)==ik)),'color',colors(ik,:)); hold on;
        end
        if slow || doSeparate,
            if ~exist('xlim','var'), xlim = get(gca,'XLim'); end
            set(gca,'XLim',xlim);
        end
        if slow,
            pause;
        end
    end
end

if exist('xlim','var'),
    set(gca,'XLim',xlim);
end


