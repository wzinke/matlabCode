function dateDiff = diffDates(in1)

flip = 0;
if ~any(size(in1) ~= 1),
    error('Not supported for matrices at this time');
elseif length(in1) == 1
    dateDiff = [];
    return
elseif size(in1,1) == 1 && size(in1,2) > 1
    in1 = in1';
    flip = 1;
end

% cell1 = cell2mat(in1(1:end),ones(1,length(
dateDiff = [];
cell1 = mat2cell(in1(1:(end-1)),ones((length(in1)-1),1));
cell2 = mat2cell(in1(2:end),ones((length(in1)-1),1));

dateDiff = cellfun(@singleDiff,cell1,cell2);

if flip,
    dateDiff = dateDiff';
end

function chk = singleDiff(day1,day2)
    chk = 0;
    if day1 > day2,
        while addDate(day2,chk) < day1
            chk = chk+1;
        end
        chk = -chk;
    elseif day2 > day1
        while addDate(day1,chk) < day2
            chk = chk+1;
        end
    end
end

end
