function klPlotSessv1(sessDate,varargin)

% Set defaults
monk = 'Gauss';
task = 'MG';
vRange = [-200:800];
mRange = [-800:200];
wvTimes = ((1:32).*25)-(9*25); 
srtHasGo = 1;
saveSess = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-m','monk'}
            monk = varargin{varStrInd(iv)+1};
        case {'-s','save'}
            saveSess = varargin{varStrInd(iv)+1};
        case {'-w','wave'}
            wave = varargin{varStrInd(iv)+1};
    end
end

% Set constants
% rootFold = 'C:\Users\Kaleb\Box Sync\Data';
rootFold = 'Y:\Users\Wolf\ephys_db';
saveDir = './Plots/sessionPlots';

% Convert date string

% Check for number of channels
chans = dir(sprintf('%s/%s/%s/DSP/DSP*',rootFold,monk,sessDate));
chanNames = {chans.name};
pullStr = @(str,inds) (str(inds));
for ic = 1:length(chanNames),
    chanStr{ic} = pullStr(chanNames{ic},[4,5]);
    chanID{ic}  = pullStr(chanNames{ic},6);
end

if max(cellfun(@str2num,chanStr)) <= 24,
    maxChan = 24;
else
    maxChan = 32;
end

maxUnit = max(cellfun(@abc2num,chanID));

% Open new figure, make axes
vSessFig = figure();
vSessAx = depthAxes('-c',maxChan,'-u',maxUnit,'-f',vSessFig);

mSessFig = figure();
mSessAx  = depthAxes('-c',maxChan,'-u',maxUnit,'-f',mSessFig);

waveFig = figure();
waveAx  = depthAxes('-c',maxChan,'-u',maxUnit,'-f',waveFig);

% Loop down axes and units to plot SDF
for ia = 1:maxChan,
    thisChan = num2str(ia); while length(thisChan) < 2, thisChan = ['0',thisChan]; end
    for iu = 1:maxUnit,
        thisUnit = lower(num2abc(iu));
        myChan = sprintf('DSP%s%s',thisChan,thisUnit);
        if ismember(myChan,chanNames),
            waveFiles = dir(sprintf('%s/%s/%s/DSP/%s/waves/*%s*.mat',rootFold,monk,sessDate,myChan,task));
            chanFiles = dir(sprintf('%s/%s/%s/DSP/%s/*%s*.mat',rootFold,monk,sessDate,myChan,task));
            if ~isempty(chanFiles),
                % Load file
                load(sprintf('%s/%s/%s/DSP/%s/%s',rootFold,monk,sessDate,myChan,chanFiles(1).name));
                
                % Get SDF
                [vSDF,vTimes] = klSpkRatev2(spiketimes,'-q',1);
                if srtHasGo,
                    [mSDF,mTimes] = klSpkRatev2(spiketimes-repmat(Task.SRT+Task.GoCue,1,size(spiketimes,2)),'-q',1);
                else
                    [mSDF,mTimes] = klSpkRatev2(spiketimes-repmat(Task.SRT,1,size(spiketimes,2)),'-q',1);
                end
                axes(vSessAx(ia,iu));
                pltMeanStd(vTimes(ismember(vTimes,vRange)),nanmean(vSDF(:,ismember(vTimes,vRange)),1),nanstd(vSDF(:,ismember(vTimes,vRange)),[],1)./sqrt(sum(isfinite(vSDF(:,ismember(vTimes,vRange))),1)),'k');
                vline(0);
                set(gca,'XTickLabel','','YTickLabel','','XLim',[vRange(1),vRange(end)]);

                axes(mSessAx(ia,iu));
                pltMeanStd(mTimes(ismember(mTimes,mRange)),nanmean(mSDF(:,ismember(mTimes,mRange)),1),nanstd(mSDF(:,ismember(mTimes,mRange)),[],1)./sqrt(sum(isfinite(mSDF(:,ismember(mTimes,mRange))),1)),'k');
                vline(0);
                set(gca,'XTickLabel','','YTickLabel','','XLim',[mRange(1),mRange(end)]);
            end
            if ~isempty(waveFiles),
                load(sprintf('%s/%s/%s/DSP/%s/waves/%s',rootFold,monk,sessDate,myChan,waveFiles(1).name));
                axes(waveAx(ia,iu));
%                 pltMeanStd(wvTimes,nanmean(wave.waves,1),nanstd(wave.waves,1)./sqrt(size(wave.waves,1)),'k');
                plot(wvTimes,nanmean(wave.waves,1),'color','k');
                set(gca,'XTickLabel','','YTickLabel','','XLim',[wvTimes(1),wvTimes(end)]);
            end
        end
    end
end

figure(vSessFig);
suptitle(sprintf('%s - %s',monk,sessDate));

figure(mSessFig);
suptitle(sprintf('%s - %s',monk,sessDate));

figure(waveFig);
suptitle(sprintf('%s - %s',monk,sessDate));

if saveSess
    saveas(vSessFig,sprintf('%s/%s/%s_StimOn.png',saveDir,monk,sessDate));
    saveas(mSessFig,sprintf('%s/%s/%s_SRT.png',saveDir,monk,sessDate));
    saveas(waveFig,sprintf('%s/%s/%s_Waves.png',saveDir,monk,sessDate));
end