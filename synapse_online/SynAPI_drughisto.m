close all; clear all; clc;

%% define channel and unit number
Channel = 9;

pico_store = 'pico';
spk_store = 'eSPK';

max_gap = 5;  % injection pulses with a spacing less than max_gap are considered as being part of the same series 

%% plot parameter
PeriInj = [-100, 200] .* 1000;  % in ms, time window around drug injection time stamp
plotWin = [ -30, 90] .* 1000;  % time window for showing plots

% PSTH parameter
binwd = 100;

%% initialize API connection
tdt = SynapseLive();
tdt.NEWONLY = 0;         % read all events in block every iteration
tdt.TIMESTAMPSONLY = 0;  % don't care what the snippets look like, just the ts

tdt.TYPE = {'snips','epocs', 'scalars'};
tdt.VERBOSE = false;

%% initialize figure
%figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1], 'Toolbar', 'none', 'Menu', 'none');
h = figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

while 1
      
    % get the most recent data
    tdt.update;
    spk = tdt.get_data(spk_store);

    pico = tdt.get_data(pico_store);
    
    InjTm = pico.onset;
    
    if(isempty(InjTm))
        continue;
    end
    
    InjTm(find(diff(InjTm) < max_gap) + 1) = [];
    
    InjTm = InjTm .* 1000;  % get it to ms
    clear tdt_pico InjState;
    
    Ntrials = length(InjTm);

    chanpos = spk.chan == Channel;

    unitlst = unique([0; double(sort(unique(spk.sortcode(chanpos))))]);
    unitlst(unitlst > 10) = [];
    Nunits = length(unitlst);
    
    %% refresh plot
    clf;
    [ha, pos] = tight_subplot(Nunits, 1, 0.01, 0.075, 0.05);

    for(u=1:length(unitlst))

        cunit = unitlst(u);

        if(cunit == 0)  % get MU
            unitpos = ismember(spk.sortcode, unitlst);
        else
            unitpos = spk.sortcode == cunit;
        end

        spkpos  = chanpos == 1 & unitpos == 1;

        % get spike times
        spks = 1000 .* spk.ts(spkpos)';

        % align data to injection times
        
        trialSPK = wz_get_SPKobj(spks, PeriInj, InjTm);

        subplot(ha(cunit+1));
        hold off;
        wz_spk_plot_hist(trialSPK, plotWin, [], [], binwd, 0);

        if(u > 1)
            ylabel('');
            %set(gca,'YTickLabel','');
        end

        if(u < Nunits)
            xlabel('');
            set(gca,'XTickLabel','');
        else
            xlabel('time after Injection [ms]');
        end
        
    end
    

    
    
    % wait until refresh
    pause(2);
    
%     % get snippet events
%     if isstruct(r)
%         if ~isnan(r.ts)
%             
%             if DO_RASTER
%                 data = TDTfilter(t.data, REF_EPOC, 'TIME', TRANGE);
%             else
%                 data = TDTfilter(t.data, REF_EPOC, 'TIME', TRANGE, 'TIMEREF', 1);
%             end
%             
%             if SORTCODE ~= 0
%                 i = find(data.snips.(EVENT).chan == CHANNEL & data.snips.(EVENT).sortcode == SORTCODE);
%             else
%                 i = find(data.snips.(EVENT).chan == CHANNEL);
%             end
%             
%             TS = data.snips.(EVENT).ts(i);
%             if isempty(TS)
%                 continue
%             end
%             
%             if DO_RASTER
%                 % match timestamp to its trial
%                 all_TS = cell(size(data.time_ranges, 2), 1);
%                 all_Y = cell(size(data.time_ranges, 2), 1);
%                 for trial = 1:size(data.time_ranges, 2)
%                     trial_TS = TS(TS >= data.time_ranges(1, trial) & TS < data.time_ranges(2, trial));
%                     all_TS{trial} = trial_TS - data.time_ranges(1, trial) + TRANGE(1);
%                     all_Y{trial} = trial * ones(numel(trial_TS), 1);
%                 end
%                 all_X = cat(1, all_TS{:});
%                 all_Y = cat(1, all_Y{:});
% 
%                 %plot raster
%                 subplot(2,1,1)
%                 hold on;
%                 plot(all_X, all_Y, '.', 'MarkerEdgeColor','k', 'MarkerSize',10)
%                 line([0 0], [1, trial-1], 'Color','r', 'LineStyle','--')
%                 axis tight;
%                 set(gca, 'XLim', [TRANGE(1), TRANGE(1)+TRANGE(2)]);
%                 ylabel('trial number')
%                 xlabel('time, s')
%                 title(sprintf('Raster ch=%d sort=%d, %d trials', CHANNEL, SORTCODE, trial))
%                 hold off;
%                 TS = all_X;
%                 subplot(2,1,2)
%             end
%             
%             % plot PSTH
%             NBINS = floor(numel(TS)/10);
%             hist(TS, NBINS);
%             hold on;
%             N = hist(TS, NBINS);
%             line([0 0], [0, max(N)*1.1], 'Color','r', 'LineStyle','--')
%             axis tight;
%             set(gca, 'XLim', [TRANGE(1), TRANGE(1)+TRANGE(2)]);
%             ylabel('number of occurrences')
%             xlabel('time, s')
%             title(sprintf('Histogram ch=%d sort=%d, %d trials', CHANNEL, SORTCODE, trial))
%             hold off;
%             drawnow
%         end
%     end
end