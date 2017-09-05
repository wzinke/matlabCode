function splWaves = klGetSplinesv1(waves)

cutWaves = klCutForSplinesv1(waves);

for iw = 1:size(cutWaves,1),
%     if isnan(cutWaves(iw,end)),
%         keyboard
%     end
    splWaves(iw,:) = spline(1:32,cutWaves(iw,:),1:.1:32);
end

keyboard