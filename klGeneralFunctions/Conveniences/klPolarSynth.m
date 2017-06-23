function synthData = klPolarSynth(r,rS,theta,thetaS,n,varargin)


% Set defaults
rType = 'gauss';
tType = 'gauss';

% Get random rs
switch rType,
    case {'g','gauss'},
        randRs = randn(n,1).*rS+r;
    case {'u','unif'},
        randRs = rand(n,1).*rS+r/2;
end

% Get random thetas
switch tType,
    case {'g','gauss'},
        randTs = randn(n,1).*thetaS+theta;
    case {'u','unif'},
        randTs = (rand(n,1)-.5).*thetaS+theta;
end


synthData(:,1) = cos(deg2rad(randTs)).*randRs;
synthData(:,2) = sin(deg2rad(randTs)).*randRs;

