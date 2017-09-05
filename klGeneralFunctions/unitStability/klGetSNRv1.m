function snr = klGetSNRv1(waves,varargin)

% Set Defaults
noiseC = 5;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'c'},
            noiseC = varargin{varStrInd(iv)+1};
    end
end

% Get average waveform and its peak-to-peak
Savg = nanmean(waves,1);
p2p = max(Savg)-min(Savg);

% Get noise type 1
residK = waves - repmat(Savg,size(waves,1),1);
noiseSpk = nanstd(residK(:));

% Noise type 2 seems unavailable from our current data...

% Get SNR: 
snr = p2p./(noiseSpk*noiseC);

%% Similarity function time
function thisSim = getSim(dist,lamb,dNorm)
%     dist = @(x,y) nansum((x-y).^2);
    thisSim = exp(-(dist.*lamb./dNorm));
end

end