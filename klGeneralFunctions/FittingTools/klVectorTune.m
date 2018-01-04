function vectAv = klVectorTune(values,locations)

[myX,myY] = pol2cart(klDeg2Rad(locations),values);
meanVals = nanmean([myX,myY],1);

% Convert back to polar
[rawTh,vectAv(2)] = cart2pol(meanVals(1),meanVals(2));
vectAv(1) = klRad2Deg(rawTh);