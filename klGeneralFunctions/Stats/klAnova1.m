function p = klAnova1(inMat,groupVar)

% Cut out nans
inMat(isnan(groupVar),:) = [];
groupVar(isnan(groupVar)) = [];

% Get means and ns
gm = nanmean(inMat,1);

uGroups = unique(groupVar);
k = length(uGroups);

nk = nan(1,k);
grpMean = nan(k,size(inMat,2));
for ik = 1:k,
	nk(ik) = sum(groupVar == uGroups(ik));
	grpMean(ik,:) = nanmean(inMat(groupVar == uGroups(ik),:),1);
end

% Now, calculate SST for each column
SST = nansum((inMat-repmat(gm,size(inMat,1),1)).^2,1);

% Now, calculate SSR
SSR = zeros(1,size(inMat,2));
for ik = 1:k,
	SSR = SSR + (nk(ik).*((grpMean(ik,:)-gm).^2));
end

% Subtract to get SSE
SSE = SST - SSR;

% Get MSR, MSE
MSR = SSR./(k-1);
MSE = SSE./(size(inMat,1)-k);

F = MSR./MSE;
% for i = 1:length(F),
%     p1(i) = 1-fcdf(F(i),(k-1),size(inMat,1)-k);
% end
p = 1-fcdf(F,ones(1,length(F)).*(k-1),ones(1,length(F)).*(size(inMat,1)-k));

