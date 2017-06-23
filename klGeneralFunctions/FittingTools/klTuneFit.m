function params = klTuneFit(sdf,locs)

visualize = 1;
tuneType  = 'gauss';

uLocs = unique(locs(isfinite(locs)));
if length(locs(isfinite(locs))) ~= length(uLocs) && size(sdf,1) == length(locs)
    maxRates = nan(size(uLocs));
    for il = 1:length(uLocs),
        maxRates(il) = nanmax(nanmean(sdf(locs == uLocs(il),:),1));
    end
else
    maxRates = sdf;
end

if length(maxRates) == length(uLocs) && size(maxRates,1) ~= size(uLocs,1),
    maxRates = maxRates';
end

adjLocs = uLocs;
midLoc = ceil((length(uLocs)+1)/2);
maxInd = find(maxRates == max(maxRates),1);
indDiff = maxInd-midLoc;
indShift = (1:length(uLocs))-indDiff;
if indDiff >= 0
    indShift(indShift < 1) = indShift(indShift < 1) + length(uLocs);
    adjLocs(1:(find(indShift == 1,1)-1)) = adjLocs(1:(find(indShift == 1,1)-1)) + 360;
elseif indDiff < 0
    indShift(indShift > length(uLocs)) = indShift(indShift > length(uLocs))-length(uLocs);
    adjLocs(find(indShift == 1,1):end) = adjLocs(find(indShift == 1,1):end) - 360;
end

shiftLocs(indShift) = adjLocs;
shiftRates(indShift) = maxRates;
shiftLocs(length(shiftLocs)+1) = shiftLocs(1)+360;
shiftRates(length(shiftRates)+1) = shiftRates(1);

uLocs(length(uLocs)+1) = uLocs(1)+360;
maxRates(length(maxRates)+1) = maxRates(1);

const = min(maxRates);
amp   = max(maxRates-const);
%normRates = (shiftRates-const)./amp;
normRates = shiftRates-const;

switch lower(tuneType)
    case {'norm','gauss','gaussian'}
        [ftMu,ftStd] = normfit(normRates);
%         ftCenter = ftMu+uLocs(maxRates == max(maxRates));
        ftCenter = ftMu+shiftLocs(shiftRates == max(shiftRates));
        fit = klMakeGauss(ftStd,'-m',ftCenter,'-x',shiftLocs(1):shiftLocs(end),'-c',const,'-a',amp);
        
        shiftX = shiftLocs(1):shiftLocs(end);
        shiftBackX = mod(shiftX,360);
        [sortX,sortInds] = sort(shiftBackX);
        shiftFit = fit(sortInds);
        
        if visualize
%             figure()
%             plot(sortX,shiftFit);
%             hold on;
% %             scatter(shiftLocs,shiftRates);
%             scatter(uLocs,maxRates);
            figure()
            ps = polar(klDeg2Rad(uLocs),maxRates);
            set(ps,'Marker','o','MarkerFaceColor','k','linestyle','none','markeredgecolor','k');
            hold on;
            polar(klDeg2Rad(sortX),shiftFit);
            keyboard
        end
        
        params.type     = 'Gaussian';
        params.mu       = ftCenter;
        params.sig      = ftStd;
        params.const    = const;
        params.amp      = amp;
        params.xvals    = shiftLocs(1):shiftLocs(end);
        params.fit      = fit;
    case {'ellipse'}
        
end
    