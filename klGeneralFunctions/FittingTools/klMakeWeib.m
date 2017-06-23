function yVals = klMakeWeib(alph,bet,varargin)

% Set defaults
gam = 1;
del = 0;
xVals = 0:100;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case '-g'
            gam = varargin{varStrInd(iv)+1};
        case '-d'
            del = varargin{varStrInd(iv)+1};
        case '-x'
            xVals = varargin{varStrInd(iv)+1};
    end
end

% Make function
weib = @(a,b,g,d,x) g-((g-d).*exp(-((x/a).^b)));

% Get values
yVals = weib(alph,bet,gam,del,xVals);