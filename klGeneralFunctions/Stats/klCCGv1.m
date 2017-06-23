function [ccg, ccgTimes, zInd] = klCCGv1(trainA, trainB, varargin)

% Set defaults
zTol = .0001;
nTens  = 0;           % How many decimal places should be considered
ccgType = 'time';
upSamp = 10;
minShift = -200;
maxShift = 200;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'shift','-s'},
            minShift = varargin{varStrInd(iv)+1}(1);
            maxShift = varargin{varStrInd(iv)+1}(2);
        case {'-t','type'},
            ccgType = varargin{varStrInd(iv)+1};
    end
end
            
% Check that trains are row vectors
if size(trainA,1) ~= 1 && any(size(trainA)) == 1,
    trainA = trainA';
end
if size(trainB,1) ~= 1 && any(size(trainB)) == 1,
    trainB = trainB';
end

switch ccgType,
    case {'time'},
        % Cut out nans
        trainA(isnan(trainA)) = [];
        trainB(isnan(trainB)) = [];

        % Make sure there aren't decimals...
        trainA  = round(trainA.*(10^nTens));
        trainB  = round(trainB.*(10^nTens));

        % Let's make a 2xt matrix of zeros where t = max spike time
        ceilA = ceil(trainA); ceilB = ceil(trainB);
        maxVal = max([ceilA,ceilB]); minVal = min([ceilA,ceilB]);

        zCorr = sum(ismember(trainA,trainB));
        for ii = 1:maxVal,
            forCorr(ii) = sum(ismember(trainA,trainB+ii));
            bacCorr(ii) = sum(ismember(trainA,trainB-ii));
        end
        ccg = [fliplr(bacCorr),zCorr,forCorr];
        ccgTimes = (-maxVal:maxVal)./(10^nTens);
        zInd = maxVal+1;
    case {'corr'}
        % Smooth inputs
        splA = spline(1:length(trainA),trainA,1:(1/upSamp):length(trainA));
        splB = spline(1:length(trainB),trainB,1:(1/upSamp):length(trainA));
        
        sVals = minShift:maxShift;
        ccgTimes = sVals;
        zInd = ceil(length(ccgTimes)/2);
        ccg = nan(1,length(sVals));
        for ib = 1:length(sVals),
            if ib < ceil(length(sVals)/2),
                thisA = splA((ceil((length(sVals)/2))-ib):end);
                thisB = splB;
                ccg(ib) = corr(thisB(1:length(thisA))',thisA');
            else
                thisA = splA;
                thisB = splB(ceil((ib-length(sVals)/2)):end);
                ccg(ib) = corr(thisA(1:length(thisB))',thisB');
            end
        end
        ccgTimes = sVals;
        zInd = ceil(length(ccgTimes)/2);
end
