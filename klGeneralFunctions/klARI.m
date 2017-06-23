function ari = klARI(in1,in2)


if nargin == 1,
    nij = in1;
else
    % Let's start by making our contingency table
    [~,~,~,~,nij] = klConting(in1,in2);
end
n = sum(nij(:));
a = sum(nij,1);
b = sum(nij,2);

term1 = 0; %sum over i and j, nij choose 2
for i = 1:numel(nij),
    if nij(i) >= 2,
        term1 = term1 + nchoosek(nij(i),2);
    end
end

% Let's do term2 (sum i ai choose 2) and term3 (sum j bj choose 2)
term2 = 0;
for i = 1:length(a),
    term2 = term2 + nchoosek(a(i),2);
end

term3 = 0;
for j = 1:length(b),
    term3 = term3 + nchoosek(b(j),2);
end

term4 = nchoosek(n,2);

% Now let's combine
numerator = term1 - ((term2*term3)/term4);
denominator = ((1/2)*(term2 + term3)) - ((term2*term3)/term4);

% And calculate ARI
ari = numerator/denominator;