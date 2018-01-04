function rgb = im2rgb(inImg,map,res,cap)

if nargin < 4
    cap = max(inImg(:));
end
if nargin < 3 || isempty(res)
    res = 256;
end


% Convert inImg to integer steps by subtracting the min, divide by range,
% then multiply by res
intImg = round(((inImg-nanmin(inImg(:)))./(cap-nanmin(inImg(:)))).*res);
intImg(intImg > res) = res;

if ischar(map)
    eval(sprintf('myMap=%s(%d);',map,res));
else
    myMap = map;
end

rgb = ind2rgb(intImg,myMap);

