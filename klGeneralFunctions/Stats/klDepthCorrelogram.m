function [outHist, outCount, xVals, randOut, randCount] = klDepthCorrelogram(varVals,chanNums,varargin)

% Set defaults
maxChan = 24;
randCorr = 1;
randType = 'shuffle';
randLoops = 100;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'rand'}
            randCorr = varargin{varStrInd(iv)+1};
        case {'rtype'}
            randType = varargin{varStrInd(iv)+1};
            if ismember(randType,{'specVals'}),
                randParams = varargin{varStrInd(iv)+2};
            end
        case {'randreps'},
            randLoops = varargin{varStrInd(iv)+1};
    end
end

% Check to make sure that varVals is numeric
switch class(varVals)
    case 'cell'
        switch class(varVals{1})
            case 'char'
                uVals = unique(varVals);
                newVarVals = nan(size(varVals));
                for iu = 1:length(uVals)
                    newVarVals(strcmp(varVals,uVals{iu})) = iu;
                end
                clear varVals; varVals = newVarVals;
            case 'double'
                varVals = cell2mat(varVals);
        end
end

% Loop through the length of varVals and see where values are the same
xVals = -(maxChan-1):(maxChan+1);
outHist = zeros(1,length(xVals));
outCount = zeros(1,length(xVals));
allCount = 0;
for i = 1:length(varVals),
    for ii = 1:length(varVals)
        if ii == i || isnan(varVals(ii)) || isnan(varVals(i)), continue; end
        thisX = chanNums(ii)-chanNums(i);
        allCount = allCount + 1;
        outCount(xVals == thisX) = outCount(xVals == thisX) + 1;
        if varVals(ii) == varVals(i),
            outHist(xVals == thisX) = outHist(xVals == thisX)+1;
        end
    end
end

if randCorr,
    for irand = 1:randLoops
        switch randType,
            case 'shuffle',
                newVals = varVals(randperm(length(varVals)));
            case 'inSess'
                newVals = varVals(randi(length(varVals),size(varVals)));
            case 'specVals'
                randVect = [];
                for iv = 1:length(randParams),
                    randVect = cat(1,randVect,ones(randParams(iv),1).*iv);
                end
                newVals = randVect(randi(length(randVect),size(varVals)));
                
        end
        [randOut{1,irand},randCount{1,irand}] = klDepthCorrelogram(newVals,chanNums,'rand',0);
    end
    
end
            