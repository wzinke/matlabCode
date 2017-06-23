function out = klCutForSplinesv1(in)

maxIn = max(in(:));
minIn = min(in(:));

out = in;

out(in==maxIn) = nan;
out(in==minIn) = nan;

out(:,1:2) = in(:,1:2);
out(:,(end-1):end) = in(:,(end-1):end);
