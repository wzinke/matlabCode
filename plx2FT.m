function [plxmat, cfg] = plx2FT(flnm)

%load first file on the list, initialize/reset config

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-code events for use by FieldTrip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % ____________________________________________________________________________ %
% %% Check if dependencies are fullfilled
% if(~exist('readPLXFileC','file'))
%     warning('readPLXFileC not found, please select now');
%     [~,PathName] = uigetfile({'*.mex*'},'readPLXFileC mex file');
%     addpath(PathName);
% end
% 
% % ____________________________________________________________________________ %
% %% check input data
% % if not specified use GUI to get the file
% if(~exist('plxfile','var') || isempty(plxfile))
%     [FileName,PathName] = uigetfile({'*.plx;*.PLX'},'Load plexon file');
%     plxfile = fullfile(PathName,FileName);
% end




fl_hdr   = ft_read_header(flnm);
fl_event = ft_read_event(flnm);
% fl_data  = ft_read_data(flnm);

cfg.dataset = flnm;
cfg.event   = fl_event;

for e=1:length(cfg.event)
    cfg.event(e).type  = 'StimOn';
end
cfg.trialdef.eventtype = 'StimOn';

%code to get events into conditionals
SearchOn = find([cfg.event.value] == 2651); %find locking event (here its search onset

% get trials for each condition, and the position of current search start (2651) in event list

TargIn_DA = SearchOn(find([cfg.event(SearchOn+32).value] == 3007 & [cfg.event(SearchOn+19).value] == 1111 & [cfg.event(SearchOn+15).value] == 4108 & [cfg.event(SearchOn+24).value] == 4300 & [cfg.event(SearchOn+22).value] == 5180 & [cfg.event(SearchOn+23).value] == 5500));
TargIn_DP = SearchOn(find([cfg.event(SearchOn+32).value] == 3007 & [cfg.event(SearchOn+19).value] == 2222 & [cfg.event(SearchOn+15).value] == 4108 & [cfg.event(SearchOn+24).value] == 4300 & [cfg.event(SearchOn+22).value] == 5180 & [cfg.event(SearchOn+23).value] == 5500));
DistIn_DA = SearchOn(find([cfg.event(SearchOn+32).value] == 3007 & [cfg.event(SearchOn+19).value] == 1111 & [cfg.event(SearchOn+15).value] == 4108 & [cfg.event(SearchOn+24).value] == 4300 & [cfg.event(SearchOn+22).value] == 5000 & [cfg.event(SearchOn+23).value] == 5680));
DistIn_DP = SearchOn(find([cfg.event(SearchOn+32).value] == 3007 & [cfg.event(SearchOn+19).value] == 1111 & [cfg.event(SearchOn+15).value] == 4108 & [cfg.event(SearchOn+24).value] == 4300 & [cfg.event(SearchOn+22).value] == 5000 & ([cfg.event(SearchOn+23).value] == 5735 | [cfg.event(SearchOn+23).value] == 5590)));
SalDistIn = SearchOn(find([cfg.event(SearchOn+32).value] == 3007 & [cfg.event(SearchOn+19).value] == 2222 & [cfg.event(SearchOn+15).value] == 4108 & [cfg.event(SearchOn+24).value] == 4300 & [cfg.event(SearchOn+22).value] == 5000 & [cfg.event(SearchOn+23).value] == 5680));

%Replace old search onset codes with condition specific codes
for t=1:length(TargIn_DA)
    cfg.event(TargIn_DA(t)).value = 7000;
end
for t=1:length(TargIn_DP)
    cfg.event(TargIn_DP(t)).value = 7001;
end
for t=1:length(DistIn_DA)
    cfg.event(DistIn_DA(t)).value = 7002;
end
for t=1:length(DistIn_DP)
    cfg.event(DistIn_DP(t)).value = 7003;
end
for t=1:length(SalDistIn)
    cfg.event(SalDistIn(t)).value = 7004;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define trials of interest and prepare data for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define trials of interest based on event codes
cfg.trialdef.eventvalue   = [7000 7001 7002 7003 7004]; % the value of the stimulus trigger
cfg.trialdef.prestim      = 0.5; % in seconds
cfg.trialdef.poststim     = 1; % in seconds
cfg = ft_definetrial(cfg);
%prepare data for relevant trials
% cfg.dftfreq   = [60-1*(1/10):(1/10):60+1*(1/10) ]; % filter out 60 hz line noise
% cfg.dftfilter = 'yes';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read in the analog data and add spike data 
plxmat = ft_preprocessing(cfg);	% call preprocessing, putting the output in 'trialdata'
% fl_spike = ft_read_spike(flnm);
% fl_data  = ft_read_data(flnm);
plxmat2 = ft_appendspike(cfg, plxmat);
% 
