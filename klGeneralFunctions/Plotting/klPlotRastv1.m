function sh = klPlotRastv1(spks,varargin)

% Set defaults
start = 1;
color = 'k';
marker = '.';
fig = gcf;
ax  = gca;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case 'color'
            color = varargin{varStrInd(iv)+1};
        case 'marker'
            marker = varargin{varStrInd(iv)+1};
        case {'figure','fig'}
            fig = varargin{varStrInd(iv)+1};
        case {'axes','ax'}
            ax = varargin{varStrInd(iv)+1};
        case 'start'
            start = varargin{varStrInd(iv)+1};
    end
end

% Make current figure/axis current
figure(fig);
axes(ax);
hold(ax,'on');

% Create y-axis matrix
yVect(1,:)  = (1:size(spks,1))+(start-1);
yVals       = repmat(yVect,size(spks,2),1);
sh          = plot(spks',yVals,'linestyle','none','color',color,'marker',marker);
