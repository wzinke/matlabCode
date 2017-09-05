function axH = depthAxes(varargin)

%  Set defaults
numChans = 24;
numUnits = 1;
width = .8;
margin = .1;
xOff = 0;
yOff = 0;

%  Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'-c','chans'}
            numChans    = varargin{varStrInd(iv)+1};
        case {'-f'}
            figNum      = varargin{varStrInd(iv)+1};
        case {'-w'}
            width       = varargin{varStrInd(iv)+1};
        case 'm'
            margin      = varargin{varStrInd(iv)+1};
        case 'x'
            xOff        = varargin{varStrInd(iv)+1};
        case 'y'
            yOff        = varargin{varStrInd(iv)+1};
        case {'-u','units'}
           numUnits = varargin{varStrInd(iv)+1};
    end
end
% Calculate the height of each axis
totSpace    = 1-(2*margin+yOff);
step        = totSpace/numChans;
wdStep      = width/numUnits;

% Open the correct figure
if ~exist('figNum','var'), 
    temp = figure(); 
    ver = version;
    verYear = str2double(ver((end-5):(end-2)));
    if verYear < 2014,
        figNum = temp;
    else
        figNum = temp.Number;
    end
    clear temp; 
end
figure(figNum);

% Loop through and create the needed axes (axis 1 is at the top and axis
% "numChans" is at the bottom
axH = nan(numChans,1);
for iu = 1:numUnits,
    for ic = 1:numChans,
        axH(ic,iu) = axes('position',[margin+xOff+wdStep*(iu-1),(1-margin)-(step*(ic-1))-(yOff),wdStep,step]);
    end
end
set(axH,'XTickLabel','','YTickLabel','');