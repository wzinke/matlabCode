function [p,r,aBar] = klRayleigh(angDist,varargin)

% Set defaults
alph = .05;
bim  = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'alpha','-a'}
            alph = varargin{varStrInd(iv)+1};
        case {'group','-g'}
            freq = varargin{varStrInd(iv)+1};
        case {'-b','bimodal','bi','bim'}
            bim = varargin{varStrInd(iv)+1};
    end
end

% Rotate angDist to be in the right orientation. Also make sure "freq" is
% in the same orientation. We're assuming that freq as an input is in the
% correct dimensions
if size(angDist,1) < size(angDist,2),
    angDist = angDist';
end
if exist('freq','var')
    if size(freq,1) < size(freq,2),
        freq = freq';
    end
end

% If we expect this to be bimodal, double angles as per Zar 26.8
if bim,
    angDist = angDist.*2;
    while any(angDist(:) >= 360),
        angDist(angDist >= 360) = angDist(angDist >= 360) - 360;
    end
end

% Calculate grouping
if exist('freq','var'),
    grpD = unique(diff(sort(angDist))); % Should only have one value!
    grpD(grpD == 0) = [];
end
% Calculate Mean r
sinAngs = sin(klDeg2Rad(angDist));
cosAngs = cos(klDeg2Rad(angDist));
if exist('freq','var'),
    sinAngs = sinAngs.*freq;
    cosAngs = cosAngs.*freq;
    mnY = nansum(sinAngs,1)./nansum(freq,1);
    mnX = nansum(cosAngs,1)./nansum(freq,1);
else
    mnY   = nanmean(sinAngs,1);
    mnX   = nanmean(cosAngs,1);
end

r   = sqrt(mnY.^2+mnX.^2);
aBar = klRad2Deg(asin(abs(mnY)./r));
if exist('grpD','var'),
    r = r.*klAdjRay(grpD);
end

if mnY < 0,
    aBar = -aBar;
end
if mnX < 0,
    aBar = 180-aBar;
end
if bim,
    aBar = aBar./2;
end
% Calculate Rayleigh's R and z
rayR = sum(~isnan(angDist),1).*r;
rayz = (rayR.^2)./sum(~isnan(angDist),1);

% Calculate Rayleigh P
rayP = @(n,R) sqrt(1+(4.*n)+(4.*((n.^2)-(R.^2))))-(1+2.*n);
p    = exp(rayP(sum(~isnan(angDist),1),rayR));