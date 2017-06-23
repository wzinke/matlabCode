% klTestAutoSortsv1

% Set artificial samples params
nSampsNz = 5000;
nSampsSig = 5000;

sigSD = 1;


% Set x dimension test range
xShift = 0:.5:4;

% Set gaussian noise, mu=0,sigma = 1
nzVals = randn(nSampsNz,2);

% Start looping through x
for ix = 1:length(xShift),
    % Print out for progress reporting...
    fprintf('Shifting X by %.1f SDs...\n\n',xShift(ix));
    
    % Set "signal" values shifted by xShift(ix)
    sigVals = (randn(nSampsSig,2).*sigSD)+repmat([xShift(ix),0],nSampsSig,1);
    
    % Sort the data
%     [sortK(ix),~,sortGap(ix,:)] = klAutoSortv1a([nzVals;sigVals]);
    
    figure(ix);
    [n,c] = hist3([nzVals;sigVals],[30,30]);
    surf(c{1},c{2},n);
%     if ix ~= length(xShift),
%         for ib = 1:length(sprintf('%.1f SDs...',xShift(ix))),
%             fprintf('\b');
%         end
%     end
end

% Because these come from 2 "distributions," plot out the second column of
% sortGap against the xShift
figure();
plot(xShift,sortGap(:,2));

% Put a vertical line in the first place that sortK = 2
vline(xShift(find(sortK==2,1)));


