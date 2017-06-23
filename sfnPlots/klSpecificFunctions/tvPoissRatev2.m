function allSpks = tvPoissRatev2(inSDF,varargin)

% Set defaults
res = 5; % In ms
nTrs = 100;
offset = -500;
visualize = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-o'},
            offset = varargin{varStrInd(iv)+1};
        case {'-r'},
            res = varargin{varStrInd(iv)+1};
        case {'-t','-n'},
            nTrs = varargin{varStrInd(iv)+1};
    end
end

% First, set up spiking function
poissRate = @(percPoint,lamb) log(1-percPoint)./(-1*lamb);

% Let's start by doing this the long way...
% Loop over inSDF by resolution to get mean firing rates
nRes =  floor(length(inSDF)./res);
myLambs = nan(1,nRes);
allITI = [];
maxSpkTm = zeros(1,nTrs);
[trSpks{1:nTrs,1}] = deal([]);
minResTms = zeros(nTrs,1);

[trSpks{1:nTrs,1}] = deal([]);
% Start trial loop
for it = 1:nTrs,
    
    cumBin = 0;
    % Start bin loop
    for ir = 1:nRes
        % Unfortunately recalculated every loop, but shouldn't be it
        % dependent
        myTimes = (1:res)+(res*(ir-1));
        myLambs(ir) = nanmean(inSDF(myTimes));
        minT = min(myTimes);
        maxT = max(myTimes);
        % Don't add spikes if nan mean SDF
        thisSectionITIs = 0;
        thisSectionTimes = 0;
        thisSectionStart = 0;
%         if minT > 525,
%             keyboard
%         end
        if ~isnan(myLambs(ir)),
            % While so we don't grab too many
            while thisSectionTimes(end) < (res),
                % Get next ITI
                randVal = rand;
                nextITI = poissRate(randVal,myLambs(ir));
                thisSectionITIs(end+1) = nextITI*1000;
                thisSectionTimes(end+1) = sum(thisSectionITIs);
%                 cumBin = cumBin + nextITI;
%                 if cumBin < res,
%                     trSpks{it} = cat(2,trSpks{it},cumBin*1000);
%                 end
            end
            while thisSectionTimes(end) > res,
                thisSectionTimes(end) = [];
            end
        end
        trSpks{it} = cat(2,trSpks{it},round(thisSectionTimes(2:end)+(res/2))+minT);
        
    end
    if mod(it,10) == 0,
        fprintf('Doing loop %d of %d...\n',it,nTrs);
    end
end

nSpks = cellfun(@length,trSpks);
allSpks = nan(nTrs,max(nSpks));
for it = 1:nTrs,
    allSpks(it,1:nSpks(it)) = round(trSpks{it});
end
allSpks = allSpks + offset;
                    
            
            
            
            
            
            