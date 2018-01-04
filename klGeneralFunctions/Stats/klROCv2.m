%% Given observed values in conditions A and B, calculate ROC curve
%  As with the poisson latency analysis, this aims to be a painstakingly
%  detailed, step-by-step process through creating ROC curves in such a way
%  that I, and hopefully others, will be able to understand the process of
%  doing this calculation.
%
%  Kaleb Lowe - 7/31/15

function [auc, roc, bootP, bootCI] = klROCv1(condA,condB,varargin)

visualize   = 0;
rocPerc     = 1:100;
allMin      = floor(min([condA(:);condB(:)]));
allMax      = ceil(max([condA(:);condB(:)]));
nSamps      = 100;
bootSig     = 1;
nBoot       = 500;
alph        = .05;

% Decode varargin here?
varStrInd = find(cellfun('isclass',varargin,'char'));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case 'boot'
            bootSig = varargin{varStrInd(iv)+1};
    end
end

res = (allMax-allMin+2)/nSamps;

if visualize
    % Start off visualizing the distributions of these two conditions
    binStep = (allMax-allMin)/100;
    bins = allMin:binStep:allMax;
    histCond(:,1) = hist(condA,bins);
    histCond(:,2) = hist(condB,bins);
    figure();
    bar(histCond);

end

% Let's cut out nans
% condA(isnan(condA)) = [];
% condB(isnan(condB)) = [];

if isempty(condA) || isempty(condB)
    auc = nan;
    roc = nan(2,nSamps);
    bootP = nan;
    bootCI = [nan,nan];
    return
end

% Let's sort the conditions in ascending order
sortA = sort(condA,'ascend');
sortB = sort(condB,'ascend');

% Conceptually, we want to make a curve for false positives vs hits. We'll
% do this by finding the value of condition B where x% of observations fall
% below that value and then calculate how many values of condition A (y%)
% also fall below that value. (I think??)
% We can implement this by sweeping from floor(allMin)-1 to ceil(allMax)+1?
rocRange = (allMin-1):res:(allMax+1-res);
for ir = 1:length(rocRange)
    percA(:,ir) = sum(sortA <= rocRange(ir))./sum(isfinite(sortA),1);
    percB(:,ir) = sum(sortB <= rocRange(ir))./sum(isfinite(sortB),1);
    if ir > 1
        % Area under segment - trapezoidal area integration between point
        % ir and ir-1. Trapezoidal area = 1/2*(h1+h2)*base
        aus(:,ir-1) = ((percB(:,ir)+percB(:,ir-1))./2).*(percA(:,ir)-percA(:,ir-1));
    end
end

% Now calculate AUC using trapezoidal integration    
auc = sum(aus,2); 
roc = [percA;percB];

% How to determine significance? For now, we'll bootstrap confidence
% intervals?
if bootSig
    % Combine all samples into one vector
    allData = [condA;condB];
    numA    = size(condA,1);
    numB    = size(condB,1);

    for ib = 1:nBoot
        % Randomize indices for allData vector
        testInd = randperm(size(allData,1));
        testA   = allData(testInd(1:numA),:);
        testB   = allData(testInd((numA+1):end),:);
        
        try
        bootAUC(:,ib) = klROCv2(testA,testB,'boot',0);      
        catch
            keyboard
        end
    end
    bootP = sum(abs(bootAUC-.5) >= abs(auc-.5),2)./nBoot;
    % Now sort the bootstrapped AUCs so we can find upper and lower tails
    % (5% total - defaults to two-tailed test)
    sortAUC = sort(bootAUC,2,'ascend');
    cutLen = floor((alph/2)*size(sortAUC,2)); % floor to be conservative
    bootCI = [sortAUC(:,cutLen),sortAUC(:,end-(cutLen-1))];
    
end


