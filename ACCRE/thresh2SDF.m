function thresh2SDF(thisFile,varargin)

rawDir = './data';
myPols = {'neg','pos'};
outDir = '/scratch/loweka/quickSpikes';

% Load behavior
load(sprintf('%s/%s/Behav.mat',rawDir,thisFile));
	
for ic = 1:64,
	
	% Load channel
	load(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',rawDir,thisFile,ic));
	
	% Figure out how big the matrix will be
	if ~isfield(Task,'trStarts'),
		Task.trStarts = trStarts;
		Task.trEnds = trEnds;
	end
	
	for ip = 1:2,
		spkTimes= chanSorts.(myPols{ip}).alltimes;
	
		nTrs= length(Task.AlignTimes);
		nSpks = nan(nTrs,1);
		for it = 1:nTrs,
			nSpks(it) = sum(spkTimes >= Task.trStarts(it) & spkTimes <= Task.trEnds(it));
		end
		maxSpk = max(nSpks(:));
		spikes.(myPols{ip}).spiketimes = nan(nTrs,maxSpk);

		% Place spikes in the matrix
		for it = 1:nTrs,
			spikes.spiketimes(it,1:nSpks(it)) =  spkTimes(spkTimes >= Task.trStarts(it) & spkTimes <= Task.trEnds(it));
		end
		spikes.(myPols{ip}).spiketimes = spikes.(myPols{ip}).spiketimes-repmat(Task.AlignTimes,1,size(spikes.(myPols{ip}).spiketimes,2));
	end
	save(sprintf('%s/%s/chan%d.mat',outDir,thisFile,ic),'spikes','-v7.3');
end

