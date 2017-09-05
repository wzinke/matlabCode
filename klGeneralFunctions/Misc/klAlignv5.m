function [outMat,currZero] = klAlignv5(inCell,indVect)


% outMat = inCell{1,1};

maxZero = max(indVect);
minZero = min(indVect);
allLengths = cellfun(@length,inCell);
maxLen = max(allLengths);
currZero = maxZero;
dim2 = maxZero+maxLen+1;
outMat = nan(size(inCell,1),dim2);
outMat2T = nan(dim2,size(inCell,1));

% Make a vector of the amount each row needs to be shifted:
shiftVect = (maxZero-indVect)+1;

% Now make this into a matrix
% Let's add the row ID...
shiftVect2 = shiftVect+([0:1:(length(shiftVect)-1)]').*dim2;

% And now increment by 1 index until the length necessary...
shiftMat = repmat(shiftVect2,1,maxLen)+repmat(0:(maxLen-1),length(shiftVect2),1);

% Make it a vector
shiftMatT = shiftMat';
shiftInd = shiftMatT(:);
% 
% % shiftMat = repmat(shiftVect+1,1,maxLen) + repmat(0:(maxLen-1),size(shiftVect,1),1);
% 
% % Now add the appropriate number to each row to make it an index
% indMat = shiftMat + repmat([0:maxLen:(maxLen*(size(shiftMat,1)-1))]',1,maxLen);
% indMatT = indMat';
% indMatVect = indMatT(:);

% Arrange the input such that it is a vector
inCellT = inCell';
inVect = cell2mat(inCellT);

% Cut out nan-correspondences from inVect, then from indMat
inVect(isnan(shiftInd)) = [];
shiftInd(isnan(shiftInd)) = [];

% Now, put that vector in the indices of outMat defined by indMat
outMat2T(shiftInd) = inVect;
outMat = outMat2T';

% % % % Loop trough the input cell
% % % for ir = 1:size(inCell,1),
% % % %     if mod(ir,1000) == 0,
% % % %         fprintf('Aligning row %d of %d\n',ir,size(inCell,1));
% % % %     end
% % %     if isnan(indVect(ir)), outMat(ir,1:size(outMat,2)) = nan(1,size(outMat,2)); continue; end
% % %     outMat(ir,(maxZero-indVect(ir)+1):(maxZero+allLengths(ir)-indVect(ir))) = inCell{ir};
% % % end
% % % while sum(isnan(outMat(:,end))) == size(outMat,1),
% % %     outMat = outMat(:,1:(end-1));
% % % end

% keyboard