function klUnitComparev1(cellA,cellB,varargin)

% Set defaults

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
	switch varargin{varStrInd(iv)},
	
	end
end

% First, get waveform cross-correlation
wvCCG = klCCGv1(cellA.wf,cellB.wf,'-t','corr');
wvSim = max(wvCCG);

% Next, get autocorrelation functions for each
autoA = klCCGv1(cellA.spikes, cellA.spikes);
autoB = klCCGv1(cellB.spikes, cellB.spikes);
autoCCG = klCCGv1(autoA,autoB,'-t','corr');
autoSim = max(autoCCG);

% Get mean rates
for ir = 1:size(cellA.spikes,1),
    rateA = sum(isfinite(cellA.spikes(ir,:))-1)./range(cellA.spikes(ir,:));
end
for ir = 1:size(cellB.spikes,1),
    rateB = sum(isfinite(cellB.spikes(ir,:))-1)./range(cellB.spikes(ir,:));
end
rateSim = abs(log(rateB)-log(rateA));