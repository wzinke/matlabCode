function [p0] = pAdjust(p)

n = length(p);
p0 = p ;

if (n <= 1)
    return
end


%%% bonferroni
%p0 = pmin(1, n * p);

%%% BH
i = n:-1:1 ;
[y,o]=sort(p,'descend');
[y,ro]=sort(o);

%baby 'cummin' function: min value up to this point in the array
for ii=1:n,
    pTemp(ii) = min(1, min(n./i(1:ii) .* p(o(1:ii))));
end
p0(1:n) = pTemp(ro);

return