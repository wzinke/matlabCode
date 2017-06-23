function [kapStar,p,kapRand] = klMatAgreev2(grp1,grp2,varargin)

% Set defaults
nRand = 10000;
doRand = 1;
randType = 'shuffle';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-n'}
            nRand = varargin{varStrInd(iv)+1};
        case {'-r'}
            doRand = varargin{varStrInd(iv)+1};
        case {'randType','randtype'},
            randType = varargin{varStrInd(iv)+1};
    end
end

% Make sure inputs are good
if length(grp1) ~= length(grp2),
    error('Unequal group lengths');
end

if size(grp1,1) ~= length(grp1),
    grp1 = grp1';
end
if size(grp2,1) ~= length(grp2),
    grp2 = grp2';
end

% Get contingency matrix
[~,~,~,exp,obs,err] = klConting(grp1,grp2);

% kapStar = sum(abs(err(:)))./length(grp1);

% Get unique labels
u1 = unique(grp1); u1(isnan(u1)) = [];
u2 = unique(grp2); u2(isnan(u2)) = [];


% Let's say kappa = sum(po-pe)/sum(1-pe) for all cells?
%   - Tested 4/19/2016, all cells = 0 (by definition possibly?)
%     Changed to same formula, but only for the largest

% Loop through
numer = 0;
denom = 0;
for i1 = 1:length(u1);
    for i2 = 1:length(u2),
        ns(i1,i2) = sum(grp1 == u1(i1) & grp2 == u2(i2));
    end
    mode1 = find(ns(i1,:) == max(ns(i1,:)),1);
    
    % Get po and pe
%     po(i1) = sum(grp1 == u1(i1) & grp2 == mode1)./length(grp1);
%     pe(i1) = (sum(grp1 == u1(i1))./length(grp1))*(sum(grp2 == mode1)./length(grp2));

    moden(i1) = ns(i1,mode1);
    pe(i1) = (sum(grp1 == u1(i1))./length(grp1))*(sum(grp2 == mode1)./length(grp2));
%     numer = numer + (po(i1) - pe(i1));
%     denom = denom + (1-pe(i1));
end

po = sum(moden)./length(grp1);
% pe = 

kapStar = (po-sum(pe))./(1-sum(pe));
kapRand = nan(1,nRand);
if doRand || nargout > 1,
    for ir = 1:nRand,
        switch randType,
            case 'shuffle'
                kapRand(ir) = klMatAgreev2(klShuffle(grp1),klShuffle(grp2),'-r',0);
            case {'random','rand'}
                kapRand(ir) = klMatAgreev2(randi(max(grp1),1,length(grp1)),randi(max(grp2),1,length(grp2)),'-r',0);
        end
    end
    p = sum(kapRand > kapStar)./length(kapRand);
end
        
        
        
        
        