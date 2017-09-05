% This function takes an input matrix (2xn or mx2) and rotates it about the
% origin by deg degrees

function out = klRotMat(in,deg)

if length(size(in)) ~= 2 || sum(size(in) == 2) < 1,
    error('Invalid input arguments');
end
turnBack = 0;
if size(in,2) > 2,
    in = in';
    turnBack = 1;
end

rad = klDeg2Rad(deg);
rotMat = [sin(rad), cos(rad); -cos(rad) sin(rad)];

out = in*rotMat;

if turnBack,
    out = out';
end