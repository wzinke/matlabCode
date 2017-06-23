function [wvMean wvStd] = getMeanWaves(waves,varargin)

% Using the input "waves", an mxn matrix with m trials and n time bins,
% calculate the mean and standard deviation at each time bin and plot
% it to the figure specified in varargin by the flag "f" and with the
% color specified by the flag "c"
%
%
%
% Created 6/18/15 by Kaleb Lowe

% Set default values
pltCol      = 'k';
doPlot      = 1;
meanType    = 'mean';
stdType     = 'sd';

% Decode varargin
if nargin > 1
    % Check that there are enough input arguments, if optional flags are
    % set
    if length(varargin{1}) ~= length(varargin), error('Not enough input arguments'); end
    for iv = 2:length(varargin{1})
        switch varargin{1}(iv)
            case 'f'
                figNum      = varargin{iv};
                figure(figNum);
            case 'c'
                pltCol      = varargin{iv};
            case 't'    
                thresh      = varargin{iv};
            case 'p'
                doPlot      = varargin{iv};
            case 'm'
                meanType    = varargin{iv};
            case 's'
                stdType     = varargin{iv};
        end
    end
end

% Get mean and standard deviation
switch lower(meanType)
    case 'mean'
        wvMean = nanmean(waves,1);
    case 'median'
        wvMean = nanmedian(waves,1);
    otherwise
        fprintf('Unrecognized summary type. Defaulting to mean\n');
end
switch stdType
    case 'sd'
        wvStd  = nanstd(waves,1);
    case 'se'
        wvStd  = nanstd(waves,1)./sqrt(size(waves,1));
    otherwise
        % Default to standard deviation
        wvStd  = nanstd(waves,1);
end

% Plot mean and standard deviation
if doPlot
    if ~exist('figNum','var') && isempty(get(0,'Children')), figure(); end
    pltMeanStd(1:size(wvMean,2),wvMean,wvStd,pltCol);
end