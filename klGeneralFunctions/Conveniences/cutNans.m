function outMat = cutNans(inMat,dim)

if nargin < 2,
    dim    = 2;
end
if dim == 1,
    inMat = inMat';
end

cut = ones(size(inMat,1),1);
for id = 1:size(inMat,2),
    cut = isnan(inMat(:,id)) & cut;
end
outMat = inMat(~cut,:);   

if dim == 1,
    outMat = outMat';
end