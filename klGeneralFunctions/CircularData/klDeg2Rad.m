function out = klDeg2Rad(in)

out = (in.*pi)./180;
out = mod(out,2*pi);