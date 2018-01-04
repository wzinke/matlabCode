function [numTradTypes,tradTypes] = klTradClassify(goodSDF,allTimeCell)

%% Get visual/movement identities (B&G-ish method)
checkVis = @(vSDF,vTimes) any(klGetConsecutive(abs(vSDF(:,vTimes >= 0 & vTimes <= 150)) > 6) >= 20);
checkMov = @(mSDF,mTimes) any(klGetConsecutive(abs(mSDF(:,mTimes <= 0 & mTimes >= -150)) > 6) >= 20);
[normVis, normMov] = klNormResp(goodSDF{1},allTimeCell{1},goodSDF{2},allTimeCell{2},'zbl');
for ii = 1:size(normVis,1)
    isVis = checkVis(normVis(ii,:),allTimeCell{1});
    [movCorr(ii),movP(ii)] = corr(normMov(ii,allTimeCell{2} >= -20 & allTimeCell{2} <= 0)',allTimeCell{2}(allTimeCell{2} >= -20 & allTimeCell{2} <= 0)');
    isMov = checkMov(normMov(ii,:),allTimeCell{2}) && (movP(ii) < .05) && (movCorr(ii) > 0);
    if isVis && isMov
        tradTypes{ii} = 'vismov';
    elseif isVis && ~isMov
        tradTypes{ii} = 'vis';
    elseif ~isVis && isMov
        tradTypes{ii} = 'mov';
    elseif ~isVis && ~isMov
        tradTypes{ii} = 'none';
    end
end
numTradTypes = nan(length(tradTypes),1);
numTradTypes(ismember(tradTypes,'vis')) = 1;
numTradTypes(ismember(tradTypes,'vismov')) = 2;
numTradTypes(ismember(tradTypes,'mov')) = 3;
numTradTypes(ismember(tradTypes,'none')) = 4;