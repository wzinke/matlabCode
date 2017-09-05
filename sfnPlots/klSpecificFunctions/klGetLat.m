function lat = klGetLat(rates,varargin)

% Set defaults
groupRates = 0;
latType    = 'avthresh';
nSig       = [2,3];
alph       = .05;
visualize  = 0;
sampRate   = 1;
raw        = 0;
base       = [-500 0];
tMax       = inf;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case '-b'
            base        = varargin{varStrInd(iv)+1};
        case '-f'
            sampRate    = varargin{varStrInd(iv)+1};
        case '-t'
            rateTimes   = varargin{varStrInd(iv)+1};
        case '-g'
            groups      = varargin{varStrInd(iv)+1};
            if length(unique(groups(isfinite(groups)))) > 1, groupRates = 1; end
        case '-type'
            latType     = varargin{varStrInd(iv)+1};
        case '-filter'
            filt        = varargin{varStrInd(iv)+1};
        case '-v'
            visualize   = varargin{varStrInd(iv)+1};
        case '-r'
            raw         = varargin{varStrInd(iv)+1};
        case '-tMax'
            tMax        = varargin{varStrInd(iv)+1};
    end
end

if raw,
    [rates,rateTimes] = klSpkRatev2(rates,'-q',1);
end

if ~exist('rateTimes','var'),
    rateTimes = 1:size(rates,2);
end

if length(rateTimes) == 2,
    rateTimes = rateTimes(1):rateTimes(2);
end

if exist('filt','var'),
    rates = rates(filt);
    if exist('groups','var'),
        groups = groups(filt,:);
    end
end

if exist('base','var'),
    blInds   = rateTimes >= min(base) & rateTimes < max(base);
    blMean   = nanmean(rates(:,blInds),2);
    blStd    = nanstd(rates(:,blInds),[],2);
    lastInd  = find(blInds,1,'last');
end

if groupRates,
    allGrps = unique(groups(isfinite(groups)));
    p = nan(1,size(rates,2));
    for ib = 1:sampRate:size(rates,2),
        p(ib) = anovan(rates(:,ib),groups,'display','off');
    end
    pAdj = nan(1,size(rates,2));
    pAdj(~isnan(p)) = pAdjust(p(~isnan(p)));
    
    firstTime = rateTimes(find(rateTimes > 0 & pAdj < alph,1));
    lat = firstTime;
        
else
    switch latType
        case 'ttest'
            % Compute paired t-test for baseline rate and post-stimulus
            % rate
            pVect = nan(1,size(rates,2));
            pAdj  = nan(1,size(rates,2));
            tMaxInd = find(rateTimes < tMax,1,'last');
            
            for ib = (lastInd+1):sampRate:tMaxInd,
                [~,pVect(ib)] = ttest(blMean,rates(:,ib));
            end
            % Adjust for multiple comparisons
            pAdj(~isnan(pVect))  = pAdjust(pVect(~isnan(pVect)));
            firstTime = rateTimes(find(pAdj < alph,1));
            
            % Put some other qualifiers here?
            
            % Set output (for now, first significant difference)
            lat = firstTime;
            if lat < 20,
                figure();
                pltMeanStd(rateTimes,nanmean(rates,1),nanstd(rates,[],1)./sqrt(size(rates,1)),'k');
                
                keyboard
            end
        case 'avthresh'
%             sig = nanstd(blMean)./sqrt(length(blMean));
%             mu  = nanmean(blMean);
%             mnRate = nanmean(rates,1);
            sig = mad(blMean,1)./sqrt(length(blMean));
            mu  = nanmedian(blMean);
            mnRate = nanmean(rates,1);
            threshInd = find(rateTimes > 0 & mnRate > (mu + nSig(1)*sig));
            if isempty(threshInd),
                lat = nan;
                
            else
                firstTime = (threshInd(1));

                if length(nSig) > 1,
                    doubleThreshInd = find(rateTimes > 0 & mnRate > (mu + nSig(2)*sig));
                    if ~isempty(doubleThreshInd),
                        if ~any(mnRate(firstTime:doubleThreshInd(1)) < (mu + nSig(1)*sig)),
                            lat = rateTimes(firstTime);
                        else
                            % Build in loop to cycle through and find lower
                            % threshold crossings
                            lat = nan; % For now, need to think on how to implement that loop
                        end
                    end
                end
                
            end
            if lat == 1,
                figure();
                pltMeanStd(rateTimes,nanmean(rates,1),nanstd(rates,[],1)./sqrt(size(rates,1)),'k');
                
                keyboard
            end
%             keyboard
        case 'trthresh'
            blThresh = blMean + nSig(1).*blStd;
            trOverThresh = nan(size(blThresh));
            trLats       = nan(size(blThresh));
            for it = (lastInd+1):sampRate:size(rates,2),
                % Check if this time index is over the threshold
                threshCheck  = rates(:,it) >= blThresh;
                trOverThresh(threshCheck == 1) = 1;
                trLats(isnan(trLats) & threshCheck == 1) = ones(sum(isnan(trLats) & threshCheck == 1),1).*rateTimes(it);
                if sum(~isnan(trOverThresh)) == size(trOverThresh,1),
                    break
                end
            end
            lat = nanmedian(trLats);
            sigTimes = lat;
            
%             keyboard
        case 'poisson'
            % Raw must be set to zero, and spiketimes must be given!
            poissLat_Longv4(rates);
        otherwise
            % Default to ttest
            % Compute paired t-test for baseline rate and post-stimulus
            % rate
            pVect = nan(1,size(rates,2));
            pAdj  = nan(1,size(rates,2));
            for ib = (lastInd+1):sampRate:size(rates,2),
                [~,pVect(ib)] = ttest(blMean,rates(:,ib));
            end
            % Adjust for multiple comparisons
            pAdj(~isnan(pVect))  = pAdjust(pVect(~isnan(pVect)));
            firstTime = rateTimes(find(pAdj < alph,1));
            
            % Put some other qualifiers here?
            
            % Set output (for now, first significant difference)
            lat = firstTime;
            %keyboard
    end
end

if visualize
    if groupRates
        colors = 'rgbcmkrgbcmkrgbcmkrgbcmkrgbcmk';
        figure(); hold on;
        for ig = 1:length(allGrps),
            pltMeanStd(rateTimes,nanmean(rates(groups == allGrps(ig),:),1),nanstd(rates(groups == allGrps(ig),:),1)./sqrt(sum(groups == allGrps(ig))),colors(ig));
        end
        vline(lat);
        
    else
        figure();
        pltMeanStd(rateTimes,nanmean(rates,1),nanstd(rates,1)./sqrt(size(rates,1)),'k');
        vline(lat);
        text(lat,1.1*max(nanmean(rates,1)+(nanstd(rates,1)./sqrt(size(rates,1)))),sprintf('Latency = %d',lat));
        
    end
end

% keyboard
            
   