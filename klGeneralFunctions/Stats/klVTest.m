function [p,h] = klVTest(angDist,tstAng)

% Set defaults
alph = .05;

if size(angDist,1) < size(angDist,2),
    angDist = angDist';
end
if size(tstAng,1) > size(tstAng,2),
    tstAng = tstAng';
end
if length(tstAng) < size(angDist,2),
    tstAng = repmat(tstAng,1,size(angDist,2));
end
% Calculate Mean r
sinAngs = sin(klDeg2Rad(angDist));
cosAngs = cos(klDeg2Rad(angDist));
mnY   = nanmean(sinAngs,1);
mnX   = nanmean(cosAngs,1);

r   = sqrt(mnY.^2+mnX.^2);
aBar = klRad2Deg(asin(abs(mnY)./r));

aBar(mnY < 0) = -aBar(mnY < 0);
aBar(mnX < 0) = 180-aBar(mnX < 0);

% Calculate Rayleigh's R and z
rayR = sum(~isnan(angDist),1).*r;
rayz = (rayR.^2)./sum(~isnan(angDist),1);

% Calculate V
V = rayR.*cos(klDeg2Rad(aBar-tstAng));
u = V.*(sqrt(2./(sum(~isnan(angDist),1))));

% u(alpha,n) approaches a one-tailed normal deviate...
p = 1-normcdf(u);
h = p < alph;