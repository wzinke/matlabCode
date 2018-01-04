function outMap = klCMapv2(pre,post,varargin)

% Assumes the inputs are 3-element color values, and the last is the #
% steps
points(1,:) = pre;
points(2,:) = post;

myIns = varargin;
while length(myIns{1}) == 3,
    points(size(points,1)+1,:) = myIns{1};
    myIns = myIns(2:end);
end
steps = myIns{1};

% Make 100 steps for each change in points
bigMap = nan((size(points,1)-1)*5000,3);
for ip = 1:(size(points,1)-1),
    range = points(ip+1,:)-points(ip,:);
    delt = range./(5000);

    % Get each element
    for ic = 1:size(pre,2),
        if delt(ic) == 0,
            bigMap((((ip-1)*5000)+1):(ip*5000),ic) = points(ip,ic);
        else
            bigMap((((ip-1)*5000)+1):(ip*5000),ic) = points(ip,ic):delt(ic):(points(ip+1,ic)-delt(ic));
        end
    end
end

numPoints = size(bigMap,1);
outStep = ceil(numPoints/(steps-1));
outMap = bigMap(1:outStep:end,:);
outMap(size(outMap,1)+1,:) = points(end,:);