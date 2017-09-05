function [p, chiSqr, df, expFreq, inMat, diffFreq, u1Save, u2Save] = klConting(in1,in2)

if isempty(in1) || isempty(in2),
    p = [];
    chiSqr = [];
    df = [];
    expFreq = [];
    inMat = [];
    diffFreq = [];
    return
end
if nargin < 2,
    inMat = in1;
else
    u1 = unique(in1);
    u1Save = u1;
    if iscell(u1(1)),
        if ischar(u1{1}),
            u1(strcmp(u1,'')) = [];
            tmp1 = nan(size(in1));
            for i1 = 1:length(u1),
                tmp1(strcmp(in1,u1{i1})) = i1;
            end
            in1 = tmp1;
            u1 = unique(tmp1);
        elseif isnumeric(u1{1}),
            u1 = cell2mat(u1);                   
            u1(isnan(u1)) = [];
        end
    end
    u1(isnan(u1)) = [];
    
    u2 = unique(in2);
    u2Save = u2;
    if iscell(u2(1)),
        if ischar(u2{1}),
            u2(strcmp(u2,'')) = [];
            tmp2 = nan(size(in2));
            for i2 = 1:length(u2),
                tmp2(strcmp(in2,u2{i2})) = i2;
            end
            in2 = tmp2;
            u2 = unique(tmp2);
        elseif isnumeric(u2{1}),
            u2 = cell2mat(u2);                   
            u2(isnan(u2)) = [];
        end
    else
        u2(isnan(u2)) = [];
    end
    
    for ir = 1:length(u1),
        for ic = 1:length(u2),
            inMat(ir,ic) = sum(in1 == u1(ir) & in2 == u2(ic));
        end
    end
end

% Get number of observations
n = sum(inMat(:));

% Calculate expected frequencies
for ir = 1:size(inMat,1),
    for ic = 1:size(inMat,2),
        expFreq(ir,ic) = (sum(inMat(:,ic))/n)*(sum(inMat(ir,:)));
    end
end
diffFreq = inMat-expFreq;

% Calculate chi squared
chiSqr  = nansum(nansum(((inMat-expFreq).^2)./expFreq));
df      = (size(inMat,1)-1)*(size(inMat,2)-1);
p       = 1-chi2cdf(chiSqr,df);