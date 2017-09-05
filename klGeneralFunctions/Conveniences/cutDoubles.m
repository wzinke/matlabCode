function out = cutDoubles(in),

% Make sure "in" is oriented correctly
dims = size(in);
if dims(1) > 1 && dims(2) == 1,
    in = in';
    flip = 1;
else
    flip = 0;
end

inMod = in;

% Shift the matrix by one
inShft = [in(2:end),nan];

out=in;
for i = 1:length(in),
    if inMod(i) == inShft(i),
        out(i+1) = nan;
    end
end

out = out(~isnan(out));
if flip
    out = out';
end
