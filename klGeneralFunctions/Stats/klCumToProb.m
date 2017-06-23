function out = klCumToProb(in)

shift = [0,in(1:(end-1))];
out = in-shift;
