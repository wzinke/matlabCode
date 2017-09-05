function [mu,std,amp] = klGaussFit(data)

% To start, let's plot the data and a first estimate
figure();
hold on;
plot(data);

% Now the initial estimate:
t=1:length(data);
halfPoint1 = find(data > (max(data)/2),1);
halfPoint2 = find(data < (max(data)/2),1,'last');
stdEst = (t(halfPoint2)-t(halfPoint1))/2;
muEst  = t(find(data == max(data),1));
ampEst = max(data);
constEst = min(data);

firstPass = klMakeGauss(stdEst,'-m',muEst,'-a',ampEst,'-t',t,'-c',constEst);
resid = nansum(firstPass-data);
plot(firstPass);

% Pick a starting direction


keyboard