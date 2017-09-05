%%normData takes a matrix "data" and normalizes it along the
%specified direction "direction." Data is normalized according
%to the input argument "normType" where normType == 1 for
%subtracting averages, and normType == 2 to compute z-scores.
%"Direction" and "normType" are optional input arguments;
%default values are 2 for direction and 1 for normType.
%normType == 3 when combined with an optional "-b" flag with
%baseline values allows for % baseline computation
%
%Currently only supports direction "2"

function adjData=normData(data,varargin)

normType = 1;
direction = 2;

if length(size(data)) > 2
    error('Currently only 2 dimensional arrays are supported');
end

if any(strncmp(varargin,'-n',2))
    normType = varargin{find(strncmp(varargin,'-n',2))+1};
end

if any(strncmp(varargin,'-d',2))
    argInd = strncmp(varargin,'-d',2);
    if length(argInd) > 1
        warning('more than one direction selected. Defaulting to first input');
        argInd = argInd(1);
    end
    direction = str2double(varargin{argInd}(find(ismember(varargin{argInd},'d'))+1));
end

if any(strncmp(varargin,'-b',2))
    argInd = find(strncmp(varargin,'-b',2));
    if length(argInd) > 1
        warning('More than one baseline vector selected. Defaulting to first input');
        argInd = argInd(1);
    end
    baselineVal = varargin{argInd+1};
end

if any(strncmp(varargin,'-z',2))
    argInd = find(strncmp(varargin,'-z',2));
    normType = 2;
    if length(argInd) > 1
        warning('More than one set of z data. Defaulting to first input');
        argInd = argInd(1);
    end
    baselineMean = varargin{argInd+1}(1);
    baselineStd  = varargin{argInd+1}(2);
end

if normType == 3 && ~any(strncmp(varargin,'-b',2))
    warning('No baseline values selected for Percent baseline. Defaulting to Z Score');
    normType = 2;
end

if normType == 2 && ~any(strncmp(varargin,'-z',2))
    warning('No Z-score mean or SD provided. Gathering from input data');
    baselineMean = nanmean(data(1:numel(data)));
    baselineStd = nanstd(data(1:numel(data)),[],1);
end

% allData = [];
% for ir = 1:size(data,1)
%     allData = cat(2,allData,data(ir,:));
% end

dataMean = nanmean(data(1:numel(data)));%nanmean(allData,2);
dataStd = nanstd(data(1:numel(data)),[],2);%allData,[],2);

switch normType
    case 0 
        adjData = data;  % this part allows for the function to be called in a script without an if statement when no adjustment is needed 
    case 1
        adjData = data - dataMean;
    case 2
        adjData = (data - baselineMean)./baselineStd;
    case 3
        adjData = (data./baselineVal).*100;
    case 4
        % Get max, min, range
        maxData = max(data(:));
        minData = min(data(:));
        dataRange = maxData-minData;
        
        % Adjustment sets minData to 0 and maxData to 1 by calculating
        % (dataPoint-minData)/range
        % If the data range is 0 (i.e., all values are the same), set all
        % values to .5
        if dataRange == 0,
            adjData = ones(size(data)).*.5;
        else
            adjData = (data-minData)./dataRange;
        end
    otherwise
        error('Unsupported normalization type');
end

% 
% 
% avData = nanmean(data,direction);
% adjData = data-repmat(avData,1,size(data,direction));
% 
% if normType == 2
%     stdData = nanstd(data,[],direction);
%     adjData = adjData./(repmat(stdData,1,size(data,direction)));
% end
% 
% 
