function wTrace = klTuneTime(sdf,locs,varargin)

% Set defaults
smoothTime = 20;
step       = 5;
times      = 1:size(sdf,2);

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)},
        case 'smooth'
            smoothTime = varargin{varStrInd(iv)+1};
        case {'res','step'}
            step = varargin{varStrInd(iv)+1};
        case {'-t','time'}
            times = varargin{varStrInd(iv)+1};
    end
end

% Loop through time

for it = 1:step:(size(sdf,2)-(smoothTime-1)),
    params = klTuneEllipse(sdf(:,(it:it+smoothTime)),locs);
    wTrace(it) = params.sig;
    aTrace(it) = params.amp;
    pdTrace(it) = params.mu;
    if wTrace(it) > 1 || wTrace(it) < 0
        wTrace(it) = nan;
    end
    if mod(it,100) == 0,
        if exist('backLen','var'), for ib = 1:backLen, fprintf('\b'); end; end
        fprintf('Time %d of %d... (%.3f%%)',it,(size(sdf,2)-(smoothTime-1)),it*100/(size(sdf,2)-smoothTime-1));
        backLen = length(sprintf('Time %d of %d... (%.3f%%)',it,(size(sdf,2)-(smoothTime-1)),it*100/(size(sdf,2)-smoothTime-1)));
    end
end

keyboard