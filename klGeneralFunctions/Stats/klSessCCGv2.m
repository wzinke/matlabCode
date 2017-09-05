function [cellCCG, times, shiftCorr] = klSessCCGv2(cellA,cellB)

% Set defaults
nTrials = 1;

% Decode varargin

% Get the raw CCG
[cellCCG, times] = ccgLoop(cellA,cellB);

% Get shift predictor by shifting B nTrials trials, then recalculating CCG
shiftB = shiftSpikes(cellB,nTrials);

% Recalculate CCG
shiftCCG = ccgLoop(cellA,shiftB);

shiftCorr = cellCCG - repmat(nanmean(shiftCCG,1),size(cellCCG,1),1);

    function [allCCG, allTimes] = ccgLoop(spikesA,spikesB)
        allCCG = [];
        for ir = 1:size(spikesA,1)
            [ccg,ccgTimes,zInd] = klCCGv1(spikesA(ir,:),spikesB(ir,:));
            if isempty(allCCG),
                allCCG = ccg; currZ = zInd;
            else
                if currZ < zInd,
                    allCCG = cat(2,nan(size(allCCG,1),zInd-currZ),allCCG,nan(size(allCCG,1),zInd-currZ));
                    currZ = zInd;
                else
                    ccg = cat(2,nan(1,currZ-zInd),ccg,nan(1,currZ-zInd));
                end
        %         keyboard
                allCCG = cat(1,allCCG,ccg);
            end
        end
        allTimes = (-currZ+1):(currZ-1);
    end

    function spikeShift = shiftSpikes(spikes,nTrials)
        spikeShift = nan(size(spikes));
        spikeShift(1:(size(spikes,1)-nTrials),:) = spikes((nTrials+1):end,:);
        spikeShift((size(spikes,1)-nTrials+1):end,:) = spikes(1:nTrials,:);
    end
end