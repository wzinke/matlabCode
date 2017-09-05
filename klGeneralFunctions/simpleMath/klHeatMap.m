function fNum = klHeatMap(heatVals,varargin)

% Set defaults
doDiag = 0;
rowName = 'x';
colName = 'y';
zName = 'z';
cmap = 'jet';
hmAx = [.1 .45 .8 .5];
diagAx   = {[.1 .1 .35 .25],[.55 .1 .35 .25]};

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-f','fig'},
            fNum = varargin{varStrInd(iv)+1};
        case {'-d','diag'},
            doDiag = varargin{varStrInd(iv)+1};
        case {'-r','rname','row'},
            rowName = varargin{varStrInd(iv)+1};
        case {'-c','cname','col'},
            colName = varargin{varStrInd(iv)+1};
        case {'-z','zname'},
            zName = varargin{varStrInd(iv)+1};
        case {'cmap'},
            cMap = varargin{varStrInd(iv)+1};
    end
end

% Open figure
if exist('fNum','var'),
    figure(fNum);
else
    fNum = figure();
    if strcmp(class(fNum),'matlab.ui.Figure'), fNum = fNum.Number; end
end

% Plot the surface
if doDiag, axes('position',hmAx); end
surf(heatVals); colormap(cMap);
ylabel(sprintf('%s',rowName)); 
xlabel(sprintf('%s',colName)); 
zlabel(sprintf('%s',zName));

% Plot the diagonal
if doDiag,
    axes('position',diagAx{1}); 
    plot(1:min(size(heatVals)),diag(heatVals));
    xlabel(sprintf('%s',rowName)); ylabel(sprintf('%s',zName));
    
    axes('position',diagAx{2});
    plot(2:min(size(heatVals)),diff(diag(heatVals(:,:))));
    xlabel(sprintf('%s',rowName)); ylabel(sprintf('Delta %s',zName));
    
end

    