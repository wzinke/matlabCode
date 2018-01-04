function [outVect, fanoVect,cvVect, skewVect,fVect]  = klCheckQuals(sdfs,times,ids)

catResp = [sdfs{1}(:,ismember(times{1},[-200:300])),sdfs{2}(:,ismember(times{2},[-300:200]))];

uIDs = nunique(ids); uIDs(uIDs==0) = [];
for ic = 1:length(uIDs)
    mySDFs = catResp(ids == uIDs(ic),:);
    grpMean = nanmean(mySDFs,1);
    outVect(ic,:) = nanvar(mySDFs,[],1)./nanvar(grpMean);
    fanoVect(ic,:) = nanvar(mySDFs,[],1)./grpMean;
    cvVect(ic,:) = nanstd(mySDFs,[],1)./grpMean;
    skewVect(ic,:) = skewness(mySDFs,[],1);
%     outVect2(ic,:) = (nanvar(mySDFs,[],1)./nanvar(grpMean)).*(sum(ids==uIDs(ic)));
end

if nargout > 4
    fVect = nan(1,size(catResp,2));
    for it = 1:size(catResp,2)
        [~,t] = anovan(catResp(:,it),ids,'display','off');
        fVect(it) = t{2,6};
    end
end