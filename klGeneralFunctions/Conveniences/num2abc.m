function out = num2abc(in)

alph = 'abcdefghijklmnopqrstuvwxyz';

switch class(in)
    case 'double'
        numFull = ceil(in/26)-1;
        if numFull == 0, fullStr = ''; else fullStr = alph(numFull); end
        numRem = mod(in,26); if numRem == 0, remStr = 'z'; else remStr = alph(numRem); end
        out = [fullStr,remStr];
end