function [lat, sigTimes] = klGetLatencyv4(spks,Task,group,varargin)

% Set defaults
alignEvents      = {'StimOnset'};
bin             = 20;
step            = 10;
rf              = 180;
alph            = .05;
rateType        = 'sdf';
sampRate        = 1;
tWind           = 20;
visualize       = 0;
parametric      = 0;
clip            = 0;

if nargin > 3
    % Check that there are enough input arguments, if optional flags are
    % set
    if length(varargin{1}) > length(varargin), 
        error('Not enough input arguments'); 
    elseif length(varargin{1}) < length(varargin),
        warning('Not enough flags set... Ignoring extra optional inputs');
    end
    for iv = 2:length(varargin{1})
        switch varargin{1}(iv)
            case 'r'
                rf = varargin{iv};
            case 'a'
                alignEvents = varargin{iv};
            case 's'
                sampRate = varargin{iv};
            case 'v'
                visualize = varargin{iv};
            case 't'
                tWind    = varargin{iv};
            case 'p'
                parametric = varargin{iv};
            case 'c'
                clip = varargin{iv};
        end
    end
end

startTic = tic;

for ia = 1:length(alignEvents)
    clear spkRate tRange spkCell aovGroup dispCell offCell pVals pVals tSub sigTimes
    alignTimes = Task.(alignEvents{ia});
    alignSpikes = spks - repmat(alignTimes,1,size(spks,2));
    
    if strcmp(alignEvents{ia},'StimOnset') && clip,
        alignSpikes(alignSpikes > repmat(Task.SRT,1,size(alignSpikes,2))) = nan;
    end
    
    switch  rateType
        case 'psth'
            [spkRate, tRange] = klSpkRatev2(alignSpikes,'-type','psth','-b',bin,'-s',step,'-q',1);
            while length(tRange) > size(spkRate,2), tRange(end) = []; end
            while length(tRange) < size(spkRate,2), tRange = cat(2,tRange,tRange(end)+step); end
        case 'sdf'
            [spkRate, tVals] = klSpkRatev2(alignSpikes,'-q',1);
            tRange           = tVals;%floor(tVals(1)):ceil(tVals(2));
    end
    
    %spkRate = spkRate(:,1:100);
    pTic = tic;
    spkCell = num2cell(spkRate,1);
    
    [aovGroup{1:size(spkRate,2)}] = deal(group);
    [dispCell{1:size(spkRate,2)}] = deal('display');
    [offCell{1:size(spkRate,2)}]  = deal('off');
    if parametric, 
        typeFun = @anovan; 
        pVals = cellfun(typeFun,spkCell(1:sampRate:end),aovGroup(1:sampRate:end),dispCell(1:sampRate:end),offCell(1:sampRate:end));
    else typeFun = @kruskalwallis;
        pVals = cellfun(typeFun,spkCell(1:sampRate:end),aovGroup(1:sampRate:end),offCell(1:sampRate:end));
    end
%     pValsKW = cellfun(@kruskalwallis,spkCell(1:sampRate:end),aovGroup(1:sampRate:end),offCell(1:sampRate:end));
    %     fprintf('pVals set in %s by v3',printTiming(pTic));
%     %% Loop and get pVals at each bin
%     pVals = nan(1,size(spkRate,2));
%     for ib = 1:size(spkRate,2),
%         pVals(ib) = anovan(spkRate(:,ib),{group},'display','off');
%     end
%     pVals   = pValsust(pVals);
    tSub   = tRange(1:sampRate:end);
%     sigTimes    = tSub(pVals < alph & tSub > nanmean(Task.StimOnset - alignTimes));
    sigTimes = tSub(klGetConsecutive(pVals < alph & tSub > nanmean(Task.StimOnset-alignTimes)) > (tWind/sampRate));
    
    if ~isempty(sigTimes), 
        lat(ia) = sigTimes(1); 
    else
        lat(ia) = nan; 
    end
    if visualize
        figure(); hold on;
        colors = 'rgbcmk';
        ug = unique(group(isfinite(group)));
        for ig = 1:length(ug)
            pltMeanStd(tRange,nanmean(spkRate(group == ug(ig),:),1),nanstd(spkRate(group == ug(ig),:),1)./sqrt(sum(isfinite(spkRate(group == ug(ig),:)),1)),colors(ig));
        end
        if ~isnan(lat(ia)), vline(lat(ia)); else hline(nanmean(spkRate(:))); end
        keyboard
    end
end

%fprintf('Latency calculated in %s\n',printTiming(startTic));
