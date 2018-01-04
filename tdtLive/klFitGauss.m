function [paramVals] = klFitGauss(targLocs,sdfVals)

targOld = targLocs;
sdfOld = sdfVals;

global targLocs sdfVals

targLocs = targOld;
sdfVals = sdfOld;

targLocs = targLocs.*45;

uLocs = nunique(targLocs);
resp = nan(1,length(uLocs));
for il = 1:length(uLocs)
    resp(il) = nanmean(sdfVals(targLocs==uLocs(il)));
end

% Get maximum value location
maxLoc = uLocs(resp==max(resp));
shiftAmnt = mod(maxLoc-180,360);

targLocs = mod(targLocs-shiftAmnt,360);

paramVals = fminsearch(@klGetResids,[min(resp),max(resp),uLocs(resp==max(resp)),45]);

paramVals(2) = mod(paramVals(2)+shiftAmnt,360);
end

function ssr = klGetResids(inParams)

global targLocs sdfVals

getGauss = @(phi,baseL,maxValL,muL,sigL) baseL + maxValL*exp(-1*.5*(((phi-muL)/sigL)^2));

uLocs = nunique(targLocs);
ssr = 0;
for il = 1:length(uLocs)
    % Get predicted value
    aPhi = getGauss(uLocs(il).*45,inParams(1),inParams(2),inParams(3),inParams(4));
    ssr = ssr + nansum((sdfVals(targLocs==uLocs(il))-aPhi).^2);
end

end