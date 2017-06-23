% if f(x) = 1-e^(-lambda*x)... then
% tv = time-varying
function x = tvPoissRate(inSDF,varargin)

% Set defaults
res = 5; % In ms
nTrs = 100;
offset = 500;
visualize = 0;

% First, set up spiking function
poissRate = @(y,lambda) log(1-y)./(-1*lambda);

% Let's start by doing this the long way...
% Loop over inSDF by resolution to get mean firing rates
nRes =  floor(length(inSDF)./res);
myLambs = nan(1,nRes);
allITI = [];
maxSpkTm = zeros(1,nTrs);
[trSpks{1:nTrs,1}] = deal([]);
minResTms = zeros(nTrs,1);
for ir = 1:nRes,
    clear spkITI spkMat
    % Get this resolution bin's mean rate (lambda for Poisson process)
    myTimes = (1:res)+(res*(ir-1));
    myLambs(ir) = nanmean(inSDF(myTimes));
    if ~isnan(myLambs(ir)),
        minT = min(myTimes);
        maxT = max(myTimes);
        % Grab trial info
        for it = 1:nTrs,
            % Get this ITI
            thisSpkTime = 0;
            startSpkTime = minResTms(it);
            while (thisSpkTime < res),
                randT = poissRate(rand,myLambs(ir))*1000;
                thisSpkTime = thisSpkTime + randT
                if thisSpkTime < res,
                    if isempty(trSpks{it}),
                        trSpks{it} = startSpkT;
                        maxSpkTm(it) = trSpks{it}(length(trSpks{it}));
                    else
                        trSpks{it} = cat(2,trSpks{it},trSpks{it}(length(trSpks{it}-1))+(poissRate(rand,myLambs(ir))*1000));
                        maxSpkTm(it) = trSpks{it}(length(trSpks{it}));
                    end
                end
            end
            spkITI(it,:) = (poissRate(rand(1,1000),myLambs(ir)));
            % Convert to times(*1000 to convert to milliseconds)
            spkMat(it,:) = cumsum(spkITI(it,:)).*1000;
            minResTms = maxT;
        end
        if visualize,
            figure();
            spkMat(spkMat > length(inSDF)) = nan;
            [sSDF,sTimes] = klSpkRatev2(spkMat);
            plot(sTimes,nanmean(sSDF,1)); hold on;
            hline(myLambs(ir));
            title(sprintf('%d - %d',minT,maxT));
            pause
            close all
        end
        spkITI(spkMat > res) = nan;
        % Cut it down so we can save some space
        goodCols = sum(isnan(spkITI),1) < size(spkITI,1);
        spkITI = spkITI(:,goodCols);    
        allITI = cat(2,allITI,spkITI);
    end
end

numSpks = (cellfun(@length,trSpks));
allSpks = nan(nTrs,max(numSpks));
for it = 1:size(trSpks,1),
    allSpks(it,1:numSpks(it)) = round(trSpks{it});
end

outSpks = nan(size(allITI));
for it = 1:size(allITI,1),
    nNotNan = sum(~isnan(allITI(it,:)));
    outSpks(it,1:nNotNan) = allITI(it,~isnan(allITI(it,:)));
end
outSpks = round(cumsum(outSpks,2).*1000);
outSpks(outSpks > length(inSDF)) = nan;
goodCols = sum(isnan(outSpks),1) < size(outSpks,1);
outSpks = outSpks(:,goodCols);

keyboard
        

