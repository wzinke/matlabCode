function distMat = klGetDist_newv1(obs,distType)

switch distType,
    case {'corr','correlation'},
        for ic = 1:size(obs,2),
            goodCols(ic) = sum(isfinite(obs(:,ic)))==size(obs,1);
        end
        distMat = 1-corr(obs(:,goodCols)',obs(:,goodCols)');
    case {'euc','sqeuc'},
        distMat = EuDist2(obs);
        if strcmp(distType,'sqeuc'),
            distMat = distMat.^2;
        end
end