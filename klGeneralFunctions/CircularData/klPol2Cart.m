function out = klPol2Cart(in)

% In (:,1) = angle (in radians), in (:,2) = radius

out(:,1) = in(:,2).*cos(in(:,1));
out(:,2) = in(:,2).*sin(in(:,1));