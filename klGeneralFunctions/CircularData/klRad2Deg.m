function out = klRad2Deg(in)

out = (in.*180)./pi;
out = mod(out,360);