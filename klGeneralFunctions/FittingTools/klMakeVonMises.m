function out = klMakeVonMises(k,varargin)

% Set defaults
mu  =0;
r   =(-pi):(pi/64):pi;
amp = 1;
const = 0;
scale = 1;
refArea = 59.7526;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case '-r',
            r = varargin{varStrInd(iv)+1};
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


vmFun = @(x,m) (1/(2*pi*besseli(0,k)))*exp(k*cos(x-m));
%vmVals = (1/(2*pi*besseli(0,k)))*exp(k*cos(r-mu));
vmVals = vmFun(r,mu);
vmAUC = integral(@(tmp) vmFun(tmp,mu),r(1),r(end));

if scale
    out = (vmVals.*((amp-const)/max(vmVals)))+const;
else
    out = (vmVals.*amp)+const;
end