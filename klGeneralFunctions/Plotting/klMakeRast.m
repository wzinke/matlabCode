%  Input:
%     spikes - nxmxc matrix where n=number of trials, m = maximum number of spikes on any given trial, and c = number of channels to do.
%              Nans should be used to fill out the matrix. That is, if on trial 1 spikes occurred at 200, 300, and 500 ms, and on trial 2 a spike occurred at 400ms,
%              this matrix would look like: [200,300,500; 400 nan nan];
%
% 	  optional inputs:
%        -f: set a figure number to plot the raster on
%        -r: a two-element vector that are the y values within which the raster will be plotted
%            that is, kMakeRast(spikes,'-r',[5,10]) will plot the resulting raster between y=5 and y=10. This is useful for superimposing a raster and SDF
%        -s: size of the plotted points

function klMakeRast(spikes,varargin),

% Set defaults
yRange = [1,size(spikes,1)];
ptSize = 1;
color = 'k';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-f','fig'},
            figNum = varargin{varStrInd(iv)+1};
        case {'-a','axes'},
            ax = varargin{varStrInd(iv)+1};
        case {'-r','range','YLim'},
            yRange = varargin{varStrInd(iv)+1};
        case {'-s','size'},
            ptSize = varargin{varStrInd(iv)+1};
        case {'-c','color'},
            color = varargin{varStrInd(iv)+1};
    end
end

% Open figure if necessary
if exist('figNum','var'), figure(figNum); end
if ~exist('figNum','var') && ~exist('ax','var'),
    figNum = figure();
end
if exist('axes','var'), axes(ax); end; hold on;

% Set up a matrix so we can scatter all at once
trCount = repmat(((1:size(spikes,1))./size(spikes,1))',1,size(spikes,2));
trCount = (trCount.*range(yRange))+yRange(1);
scatter(spikes(:),trCount(:),ptSize,color,'o','filled','markeredgecolor','none');
