function [wvWidth,tfRepol] = klWvWidthv2(alWaves,alTimes)

visualize = 0;
endCutTime = 50;

% Cut down waves and times to times > 0 and below end-endCutTime
timeCrit = alTimes > 0 & (alTimes < (max(alTimes)-endCutTime));
checkWaves = (alWaves(:,timeCrit)-repmat(nanmean(alWaves,2),1,sum(timeCrit)))./repmat(nanstd(alWaves,[],2),1,sum(timeCrit));
checkTimes = alTimes(:,timeCrit);

% Loop through waveforms
for iw = 1:size(alWaves,1),
    
    if sum(isnan(checkWaves(iw,:))) == size(checkWaves,2),
        wvWidth(iw) = nan; tfRepol(iw) = nan;
        continue
    end
    
    % Get extrema
    [xmax,imax,xmin,imin] = extrema(checkWaves(iw,:));
    imax(ismember(imax,[1,size(checkWaves,2)])) = [];
    imin(ismember(imin,[1,size(checkWaves,2)])) = [];
    
    % Get the next extrema, dependent on whether this is a negative or
    % positive spike
    if (checkWaves(iw,1) - checkWaves(iw,10)) < 0,
        if isempty(imax), wvWidth(iw) = nan; tfRepol(iw) = nan; continue; end
        
        [sortMaxInds, sortInds] = sort(imax);
        sortMaxVals = xmax(sortInds);
        wvWidth(iw) = checkTimes(sortMaxInds(sortMaxVals==max(sortMaxVals)));
        peakVal = alWaves(iw,sortMaxInds(sortMaxVals==max(sortMaxVals))+sum(alTimes <= 0));
        
        % Get the next time the wave < .75*peakVal
        tfFind = find(alWaves(iw,:) < peakVal*.75 & alTimes > wvWidth(iw));
        if isempty(tfFind),
            tfRepol(iw) = nan;
        else
            tfRepol(iw) = alTimes(tfFind(1))-wvWidth(iw);
        end
    else
        if isempty(imin), wvWidth(iw) = nan; tfRepol(iw) = nan; continue; end
        [sortMinInds, sortInds] = sort(imin);
        sortMinVals = xmin(sortInds);
        wvWidth(iw) = checkTimes(sortMinInds(sortMinVals==min(sortMinVals)));
        peakVal = alWaves(iw,sortMinInds(sortMinVals==min(sortMinVals))+sum(alTimes <= 0));
        
        % Get the next time the wave < .75*peakVal
        tfFind = find(alWaves(iw,:) > peakVal*.75 & alTimes > wvWidth(iw));
        if isempty(tfFind),
            tfRepol(iw) = nan;
        else
            tfRepol(iw) = alTimes(tfFind(1))-wvWidth(iw);
        end
    end
    
    if wvWidth(iw) < .01,
%         keyboard
    end
    
    if isnan(tfRepol(iw)),
%         keyboard
%         visualize = 1;
    end
    
    if visualize,
        x=figure();
        plot(alTimes,alWaves(iw,:));
        vline(wvWidth(iw));
        if ~isnan(tfRepol(iw)), vline(tfRepol(iw)+wvWidth(iw)); end;
        
        keyboard
        
        close(x)
    end
    
%     visualize = 0;
end

    