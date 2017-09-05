function [normVis, normMov] = klNormResp(visAlign,visTimes,movAlign,movTimes,respNormType,varargin)

blWind = -300:-100;
visWind = -200:300;
movWind = -300:200;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'bl'},
            blWind = varargin{varStrInd(iv)+1};
        case {'v','vis'},
            visWind = varargin{varStrInd(iv)+1};
        case {'m','mov'},
            movWind = varargin{varStrInd(iv)+1};
    end
end


switch respNormType,
    case {'zbl','z'},
        blMeans = nanmean(visAlign(:,ismember(visTimes,blWind)),2);
        blStds = nanstd(visAlign(:,ismember(visTimes,blWind)),[],2);
        normVis = (visAlign-repmat(blMeans,1,size(visAlign,2)))./repmat(blStds,1,size(visAlign,2));
        normMov = (movAlign-repmat(blMeans,1,size(movAlign,2)))./repmat(blStds,1,size(movAlign,2));
    case {'ztr'},
        catVals = [visAlign(:,ismember(visTimes,visWind)),movAlign(:,ismember(movTimes,movWind))];
        blMeans = nanmean(catVals,2);
        blStds = nanstd(catVals,[],2);
        normVis = (visAlign-repmat(blMeans,1,size(visAlign,2)))./repmat(blStds,1,size(visAlign,2));
        normMov = (movAlign-repmat(blMeans,1,size(movAlign,2)))./repmat(blStds,1,size(movAlign,2));
    case {'ztrbl'},
        catVals = [visAlign(:,ismember(visTimes,visWind)),movAlign(:,ismember(movTimes,movWind))];
        blMeans = nanmean(visAlign(:,ismember(visTimes,blWind)),2);
        blStds = nanstd(catVals,[],2);
        normVis = (visAlign-repmat(blMeans,1,size(visAlign,2)))./repmat(blStds,1,size(visAlign,2));
        normMov = (movAlign-repmat(blMeans,1,size(movAlign,2)))./repmat(blStds,1,size(movAlign,2));
    case {'max'},
        blMeans = nanmean(visAlign(:,ismember(visTimes,blWind)),2);
        catVals = [visAlign(:,ismember(visTimes,visWind)),movAlign(:,ismember(movTimes,movWind))];
        absMax = max(abs(catVals-repmat(blMeans,1,size(catVals,2))),[],2);
        normVis = (visAlign-repmat(blMeans,1,size(visAlign,2)))./repmat(absMax,1,size(visAlign,2));
        normMov = (movAlign-repmat(blMeans,1,size(movAlign,2)))./repmat(absMax,1,size(movAlign,2));
    case {'scale'},
        catVals = [visAlign(:,ismember(visTimes,visWind)),movAlign(:,ismember(movTimes,movWind))];
        minVals = min(catVals,[],2);
        maxVals = max(catVals,[],2);
        normVis = (visAlign-repmat(minVals,1,size(visAlign,2)))./repmat(maxVals-minVals,1,size(visAlign,2));
        normMov = (movAlign-repmat(minVals,1,size(movAlign,2)))./repmat(maxVals-minVals,1,size(movAlign,2));
    case {'base','bl'},
        blMeans = nanmean(visAlign(:,ismember(visTimes,blWind)),2);
        normVis = visAlign-repmat(blMeans,1,size(visAlign,2));
        normMov = movAlign-repmat(blMeans,1,size(movAlign,2));
    case {'none'},
        normVis = visAlign;
        normMov = movAlign;
end
