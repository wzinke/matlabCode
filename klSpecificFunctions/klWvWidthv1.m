function [wvWidth, success] = klWvWidthv2(wave,varargin)
%% Get waveform width
%
% Set defaults
wvTime      = ((1:length(wave))-9).*.25;
splTimes    = ((1:.1:32).*25)-(9*25);
minWdth     = 75;
zeroTol     = .001;
visualize   = 0;
success     = 1;
upSamp      = 1;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-t','time'},
            wvTime = varargin{varStrInd(iv)+1};
        case {'-u'},
            upSamp = varargin{varStrInd(iv)+1};
        case {'-v'},
            visualize = varargin{varStrInd(iv)+1};
    end
end

% Upsample the time and the wave by 10
if upSamp,
    splWave = spline(1:32,wave,1:.1:32);
    if length(wvTime)==length(wave),
        splTimes = spline(1:32,wvTime,1:.1:32);
    end
else
    splWave = wave;
    splTimes = wvTime;
end

% Align the wave using klAlign
[alignWv, time] = klTroughAlign(splWave,splTimes,0);

% Now, clip the wave to times > minimum width
clpWv = alignWv(time > minWdth);
clpTm = time(time > minWdth);

% Get the derivative of the wave...
divWv = abs(abs(diff(clpWv)) - diff(clpWv));
if divWv(1) > 0,
    endT = clpTm(find(divWv < zeroTol,1));
else
    endT = clpTm(find(divWv > 0,1));
end
if isempty(endT),
    endT = max(clpTm);
    success = 0;
end
wvWidth = endT;

if visualize
    figure();
    plot(time,alignWv);
    vline(0); vline(endT);
    keyboard
end




