function r2 = klGetR2(aovTable)

% Find "Error" and "Total" Rows
errRow = find(strcmp(aovTable(:,1),'Error'),1);
totRow = find(strcmp(aovTable(:,1),'Total'),1);

for iv = 2:(errRow-1),
    ssb(iv-1,1) = aovTable{iv,2};
end
sst = aovTable{totRow,2};

r2 = ssb./sst;