function out = klEucDist(x,y)

if nargin < 2,
    y = x;
end

if size(x,2) ~= size(y,2)
    error('Dimension mismatch');
end

runSum = zeros(size(x,1),size(y,1));
for ic = 1:size(x,2),
    runSum = runSum + (abs(bsxfun(@minus,x(:,ic),y(:,ic)')).^2);
end

out = sqrt(runSum);