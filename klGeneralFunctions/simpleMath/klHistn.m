function outMat  = klHistn(obs,varargin)

% Set defaults
nBins = ones(1,size(obs,2)).*30;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-b','bins'},
            nBins = varargin{varStrInd(iv)+1};
            if length(nBins) ==1,
                nBins = ones(1,size(obs,2)).*nBins;
            end
    end
end

% Get range and bin locations
minObs = min(obs,[],1);
maxObs = max(obs,[],1);

% Loop through columns
outPos = nan(size(obs));
for ic = 1:size(obs,2),
    outPos(:,ic) = 1+floor((obs(:,ic)-minObs(ic))./((maxObs(ic)-minObs(ic))./(nBins(ic)-1)));
%     outPos(:,ic) = 1+ceil((obs(:,ic)-repmat(minObs(ic),size(obs,1),1))./((maxObs(ic)-minObs(ic))./(nBins(ic)-1)));
end

% Make a string we can decode
obsStr = cell(size(obs,1),1);
[obsStr{:}] = deal('');
for ir = 1:size(obs,1),
    for ic = 1:size(obs,2),
        obsStr{ir} = cat(2,obsStr{ir},sprintf('/%d',outPos(ir,ic)));
    end
end

% Get unique combos
uCombos = unique(obsStr);

outMat = zeros(nBins);
% Loop and count
for iu = 1:length(uCombos),
    
    % Decode string
    slashInd = strfind(uCombos{iu},'/');
    evalIn = 'outMat(';
    for is = 1:length(slashInd),
        if is < length(slashInd),
            evalIn = sprintf('%s%d',evalIn,str2double(uCombos{iu}((slashInd(is)+1):(slashInd(is+1)-1))));
            evalIn = [evalIn,','];
        else
            evalIn = sprintf('%s%d',evalIn,str2double(uCombos{iu}((slashInd(is)+1):end)));
            evalIn = [evalIn,') = '];
        end
    end
    evalOut = sprintf('%d;',sum(ismember(obsStr,uCombos{iu})));
%     if sum(ismember(obsStr,uCombos{iu})) > 2,
%         keyboard
%     end
    eval(sprintf('%s%s',evalIn,evalOut));
end
