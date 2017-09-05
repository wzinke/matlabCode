function params = klTuneEllipse(sdf,locs)

visualize = 0;

% Get the max rates for the different locations
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

% Convert the locations/rates, as would be displayed in polar coordinates,
% and convert to cartesian display
cartRates = klPol2Cart([klDeg2Rad(uLocs),maxRates]);

% Fit an ellipse
try
    [cartC, cartX, cartY, cartAng] = fitellipse(cartRates);
catch
    params.type     = 'Ellipse';
    params.mu       = nan;
    params.amp      = nan;
    params.major    = nan;
    params.minor    = nan;
    params.sig      = nan;
    params.fit      = nan;
    return
end
cartEll = klMakeEllipse(cartX,cartY,'-c',cartC,'-a',cartAng);

% Convert back to polar
fitEll  = klCart2Pol(cartEll);

if visualize,
    figure()
    pRates = polar(klDeg2Rad(uLocs),maxRates);
    set(pRates,'linestyle','none','marker','o','markerfacecolor','k','markeredgecolor','k');
    hold on;
    polar(fitEll(:,1),fitEll(:,2));
    keyboard
end

params.type     = 'Ellipse';
params.mu       = klRad2Deg(fitEll(fitEll(:,2) == max(fitEll(:,2)),1));
params.amp      = max(fitEll(:,2));
params.major    = max([cartX,cartY]);
params.minor    = min([cartX,cartY]);
params.sig      = params.minor/params.major;
params.fit      = fitEll;
if length(params.mu) > 1,
    params.mu = params.mu(1);
end