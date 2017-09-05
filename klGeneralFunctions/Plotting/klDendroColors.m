function klDendroColors(inH,newColors,varargin)

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-o','old'},
            oldColors = varargin{varStrInd(iv)+1};
    end
end


currColors = cell2mat(get(inH,'color'));
uColors = unique(currColors, 'rows');

if exist('oldColors','var'),
    goodRows = [];
    for i = 1:size(oldColors,1),
        goodRows = cat(2,goodRows,find(sum(repmat(oldColors(i,:),size(uColors,1),1)==uColors,2)==3));
    end
    uColors = [nan(1,3); uColors(goodRows,:)];
end

for ic = 2:size(uColors,1),
    idxInds = find(ismember(currColors,uColors(ic,:),'rows'));
    for ii = 1:length(idxInds),
        set(inH(idxInds(ii)), 'color', newColors(ic-1,:));
    end
end