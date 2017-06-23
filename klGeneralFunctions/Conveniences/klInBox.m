function outLog = klInBox(x,y,c,w)

if numel(w) == 1,
    w(2) = w(1);
end

shftX = x-c(1);
shftY = y-c(2);

inX = shftX >= -w(1)/2 & shftX <= w(1)/2;
inY = shftY >= -w(2)/2 & shftY <= w(2)/2;

outLog = inX & inY;
