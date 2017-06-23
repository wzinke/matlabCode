function [out, xvals] = klMakeGauss(sd,varargin)

% Set defaults
amp = 1;
mu = 0;
xvals = -100:100;
scale = 1;
const = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case '-x',
            xvals = varargin{varStrInd(iv)+1};
        case '-m',
            mu = varargin{varStrInd(iv)+1};
        case '-a',
            amp = varargin{varStrInd(iv)+1};
        case '-s',
            scale = varargin{varStrInd(iv)+1};
        case '-c',
            const = varargin{varStrInd(iv)+1};
    end
end

normFun = ((1/sqrt(2*pi*sd.^2))*exp(-(xvals-mu).^2/(2*sd.^2)));
if scale
    out = (normFun.*(amp/max(normFun)))+const;
else
    out = (normFun.*amp)+const;
end