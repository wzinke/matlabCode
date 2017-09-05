function [saccStart, saccEnd] = klGetSacc(x,y,varargin)


smoothWind = 20;
saccWind   = [-inf,inf];
tVals      = 1:(length(x)-smoothWind);
accSearchWind = 100;

% Start by smoothing X and Y values (saccades should survive this smoothing process)
smX = [];
smY = [];
for ii = 1:(length(x)-(smoothWind)),
    smX(ii) = nanmean(x(ii:(ii+smoothWind)));
    smY(ii) = nanmean(y(ii:(ii+smoothWind)));
end

% Get velocities - combine to get speed (unsigned)
vX = diff(smX); vY = diff(smY);
combV = sqrt(vX.^2 + vY.^2);
combAcc = diff(combV);

velT  = tVals(1:(end-1));
% Get extrema for the velocity in the desired time window
% [maxV,maxVInd,minV,minVInd] = extrema(combV(velT >= saccWind(1) & velT <= saccWind(2)));

maxV = max(combV); maxVInd = find(combV(velT >= saccWind(1) & velT <= saccWind(2)) == max(combV(velT >= saccWind(1) & velT <= saccWind(2))));
preAcc = combAcc((maxVInd-accSearchWind):maxVInd);

checkA = inf;
for iw = 1:length(preAcc),
    if abs(preAcc((end-(iw-1)))) < checkA,
        checkA = abs(preAcc((end-(iw-1))));
    else
        break
    end
end

    

keyboard
