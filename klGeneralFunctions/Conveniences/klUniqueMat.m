function uMat = klUniqueMat(inMat,varargin)

inCell = mat2cell(inMat,ones(size(inMat,1),1),size(inMat,2));
uVals = cellfun(@nunique,inCell,'UniformOutput',0);

[histVals, histC] = cellfun(@hist,inCell,uVals,'UniformOutput',0);

uCell = cellfun(@(vals,hC,hV) ismember(vals,hC(hV==1)),inCell,histC,histVals,'UniformOutput',0);

uMat = cell2mat(uCell);