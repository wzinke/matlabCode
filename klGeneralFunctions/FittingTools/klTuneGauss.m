function fit = klTuneGauss(xPts,yVals,varargin)

visualize = 0;

varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-v','vis'},
            visualize = varargin{varStrInd(iv)+1};
    end
end

% Check that dimensions are the same (if not, rotate y)
if size(xPts,1) ~= size(yVals,1)
    if size(xPts,1) == size(yVals,2),
        yVals = yVals';
    else
        error('Inconsistent Dimensions');
    end
end

% Define a gaussian function
gaussFun = @(x,mu,std,amp,base) (exp((-(x-mu).^2)./(2*std^2)).*(amp-base))+base;
gFitFun = @(params,x) gaussFun(x, params(1), params(2), params(3), params(4));

% Get rid of nans
nanVect = isnan(xPts) | isnan(yVals);
xPts = xPts(~nanVect); yVals = yVals(~nanVect);


% Estimate mu, std, amp, and base
muEst = xPts(yVals == max(yVals)); if length(muEst) > 1, muEst = muEst(1); end
ampEst = yVals(yVals == max(yVals)); if length(ampEst) > 1, ampEst = ampEst(1); end
baseEst = yVals(yVals == min(yVals)); if length(baseEst) > 1, baseEst = baseEst(1); end;
minX = xPts(yVals == min(yVals)); if length(minX) > 1, minX = minX(1); end
stdEst = (muEst-minX)./2;

startVect = [muEst,abs(stdEst),ampEst,baseEst];

fitVect = nlinfit(xPts,yVals,gFitFun,startVect);

if fitVect(1) < -360 || fitVect(1) > 360,
    fit.mu = nan; fit.sig = nan; fit.amp = nan; fit.bl = nan;
    warning('Fit falls outside the range of possible input, suggesting a poor fit...');
else
    fit.mu = fitVect(1); while fit.mu < 0, fit.mu = fit.mu + 360; end
    fit.sig = fitVect(2);
    fit.amp = fitVect(3);
    fit.bl = fitVect(4);
end

if visualize
    figure()
    plot(xPts,yVals,'r*','linestyle','none');
    hold on;
    plot(xPts(1):10:xPts(end),gaussFun(xPts(1):10:xPts(end),fitVect(1),fitVect(2),fitVect(3),fitVect(4)),'b');
    keyboard
end