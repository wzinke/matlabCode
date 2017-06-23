function out = abc2num(in)

alph = 'abcdefghijklmnopqrstuvwxyz';

switch class(in)
    case 'char'
        mod = find(ismember(alph,lower(in(end))));
        if length(in) == 2
            mult = find(ismember(alph,lower(in(1))));
        elseif length(in) == 1
            mult = 0;
        else
            error('Unsupported input length');
        end
        out = mod + mult*26;
    case 'double'
end