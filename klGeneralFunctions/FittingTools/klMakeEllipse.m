function pnts = klMakeEllipse(xR,yR,varargin)

% Set defaults    
center = [0,0];
nSteps = 100;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case '-c',
            if length(varargin{varStrInd(iv)}) > 1,
                center = varargin{varStrInd(iv)+1};
            end
        case '-a',
            rotAng = varargin{varStrInd(iv)+1};
    end
end
     
nPnts = 0;
for x = -xR:(2*xR/(nSteps-1)):xR,
    nPnts = nPnts+1;
    rpnts(nPnts,1) = x;
    rpnts(nPnts,2) = real(sqrt((yR^2)*(1-((x^2)/(xR^2)))));
%     if ~isreal(rpnts(nPnts,2)),
%         keyboard
%     end
end

for x = xR:-(2*xR/(nSteps-1)):-xR,
    nPnts = nPnts+1;
    rpnts(nPnts,1) = x;
    rpnts(nPnts,2) = real(-sqrt((yR^2)*(1-((x^2)/(xR^2)))));
end

if exist('rotAng','var'),
    rpnts = rpnts*([cos(rotAng),-sin(rotAng); sin(rotAng),cos(rotAng)]');
end

pnts(:,1) = rpnts(:,1)+center(1);
pnts(:,2) = rpnts(:,2)+center(2);

