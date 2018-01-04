function fitVals = klFitGauss(xVal,values)

indiv = 0;

% Get mean values
[uX,uV] = umeans(xVal,values);

% Rotate everything to put max(uV) at 180
rotVal = 180-uX(uV==max(uV));
xShift = mod(xVal+rotVal,360);

% Initialize variables
init(1) = randi([0,360]);
init(2) = randi([3,10]);
init(3) = min(uV);
init(4) = max(uV);

f0 = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0,0],'Upper',[inf,inf,360,inf]);
ft = fittype('a+b*exp(-((x-c)/d)^2)','options',f0,'independent',{'x'});
f = fit(xShift,values,ft);

fitVals = coeffvalues(f);
% opts = optimset('Display','iter');
% 
% if indiv
%     fit = fminsearch(@(x) klDoFit(x,xShift,values),init,opts);
% else
%     fit = fminsearch(@(x) klDoFit(x,mod(uX+rotVal,360),uV),init,opts);
% end


end


function sse = klDoFit(params,xVal,values)

% Predict values
pred = params(3) + params(4).*exp(-.5.*(((xVal-params(1))./params(2)).^2));

% Get residuals
resid = (values-pred).^2;
sse = nansum(resid);

end

function [uVals,uMeans] = umeans(val,obs)

uVals = nunique(val);
uMeans = nan(size(uVals));
for iv = 1:length(uVals)
    uMeans(iv) = nanmean(obs(val==uVals(iv)));
end
end