function [p,x2] = klCompareX2(obs1,obs2,id1,id2,varargin)

% Set defaults
alph = .05;
xVal1 = 1:size(obs1,2);
xVal2 = 1:size(obs2,2);

% Decode varargin

% Get unique IDs
u1 = unique(id1); u1(isnan(u1)) = [];
u2 = unique(id1); u2(isnan(u2)) = [];

% Create contingency table, run statistic
for i1 = 1:length(u1),
    for i2 = 1:length(u2),
        obsMat(i1,i2) = sum(id1 == u1(i1) & id2 == u2(i2));
    end
end
[p,x2,expMat] = klConting(obsMat);

if p < alph,
    mnX2 = x2./numel(obsMat);
    sigMat = (((obsMat-expMat).^2)./expMat) > mnX2;
    for i1 = 1:length(u1), 
        for i2 = 1:length(u2),
            if sigMat(i1,i2),
                figure();
                subplot(1,2,1); pltMeanStd(xVal1,nanmean(obs1(id1 == i1,:),1),nanstd(obs1(id1 == i1,:),[],1)./sum(id1 == i1),'k');
                subplot(1,2,2); pltMeanStd(xVal2,nanmean(obs2(id2 == i2,:),1),nanstd(obs2(id2 == i2,:),1)./sum(id2 == i2),'k');
                if obsMat(i1,i2) > expMat(i1,i2), suptitle('More'); else suptitle('Fewer'); end
            end
        end
    end
end