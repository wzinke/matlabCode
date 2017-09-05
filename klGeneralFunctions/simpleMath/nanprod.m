function out = nanprod(in,dim)

% Currently only supports dim <= 2
if nargin < 2,
    dim = 1;
end

if isempty(dim) || dim > 2,
    error('Invalid dim argument');
end

flip = 0;
if dim ==2,
    flip = 1;
    in = in';
end

out = nan(size(in,1),1);
for ix = 1:size(in,1),
    out(ix) = prod(in(ix,~isnan(in(ix,:))));
end

if flip, out = out'; end;