function distMat = makeDistMat(obs,type,catVars,print)
    if nargin < 4,
        print = 0;
    end
    % Rescale again for product distance
    if ismember(type,{'prod'}),
        maxObs = nanmax(obs,[],1);
        minObs = nanmin(obs,[],1);
        normObs = (obs-repmat(minObs,size(obs,1),1))./repmat((maxObs-minObs),size(obs,1),1);
        normRaw = normObs;
    end

    % Start by computing pair-wise distances
    if print,
        fprintf('Initializing pairwise distances...')
    end
    if ismember(type,{'corr'}),
        % Clean up obs, if "corr"
        goodCols = zeros(1,size(obs,2));
        for ic = 1:size(obs,2),
            goodCols(ic) = sum(~isnan(obs(:,ic))) == size(obs,1);
        end
        obs = obs(:,logical(goodCols));
        distMat = 1-corr(obs',obs');
        for ir = 1:size(obs,1),
            for ic = ir:size(obs,1),
                distMat(ic,ir) = nan;
            end
        end
    elseif ismember(type,{'euc'}),
        distMat = EuDist2(obs(:,~catVars));
        for ir = 1:size(obs,1),
            for ic = ir:size(obs,1),
                distMat(ic,ir) = nan;
            end
        end
    else
        distMat = nan(size(obs,1),size(obs,1));
        for im = 1:size(obs,1),
            for in = (im+1):size(obs,1),
                if ismember(type,{'prod','exprod'}) || doNorm,
                    distMat(im,in) = klGetDistv2(normObs(im,:),normObs(in,:),'-t',type,'cat',catVars);
                else
                    distMat(im,in) = klGetDistv2(obs(im,:),obs(in,:),'-t',type,'cat',catVars);
                end
            end
        end
    end
    if print,
        fprintf('Done!\n');
    end
end
