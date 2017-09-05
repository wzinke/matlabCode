function outDist = klGetDistv2(obs1,obs2,varargin)

% Set defaults
type = 'euc';
catVars = zeros(1,length(obs1));
weights = ones(1,length(obs1));

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-t','type'}
            type = varargin{varStrInd(iv)+1};
        case {'cat'}
            catVars = varargin{varStrInd(iv)+1};
        case {'-w'},
            weights = varargin{varStrInd(iv)+1};
    end
end

switch type,
    case 'corr',
        nanCrit = ~isnan(obs1) & ~isnan(obs2);
        if sum(nanCrit) == 0,
            outDist = nan;
        else
            outDist = 1-corr(obs1(nanCrit)',obs2(nanCrit)');
        end
    case 'euc'
%             outDist = sqrt(nansum((obs1-obs2).^2));
        allSim = nan(1,size(obs1,2));
        if sum(catVars) > 0,
            allSim(catVars) = abs(double(obs1(catVars)) - double(obs2(catVars))) < .00001;
        end
        allSim(~catVars) = abs(obs1(~catVars) - obs2(~catVars));
        allSim = allSim.*weights;
        outDist = sqrt(nansum((allSim(~catVars).^2)));
%         if sum(catVars) > 0, outDist = outDist/nanprod(allSim(catVars)); end
        if sum(catVars) > 0, outDist = outDist.*nanprod(allSim(catVars)); end
    case 'sqeuc',
        outDist = (nansum(((obs1-obs2).^2).*weights));
    case 'prod'
        allSim = nan(1,size(obs1,2));
        if sum(catVars) > 0,
            allSim(catVars) = abs(double(obs1(catVars)) - double(obs2(catVars))) < .00001;
        end
        allSim(~catVars) = (1-abs(obs1(~catVars) - obs2(~catVars))).*weights(~catVars);
%                     thisDist = 1-nanprod(1-allDist);
%         allSim = allSim.*weights;
        outDist = 1-(nanprod(allSim)^(1/length(allSim)));
    case 'exprod',
        allSim = nan(1,size(obs1,2));
        if sum(catVars) > 0,
            allSim(catVars) = obs1(catVars) ~= obs2(catVars);
        end
        allSim(~catVars) = (abs(obs1(~catVars) - obs2(~catVars))).*weights(~catVars);
        outDist = 1-nanprod(exp(-allSim));
end
end