% Make anti events

function Task = klGetAntiTaskv2(inEvs,inEvTms,events,eyeX,eyeY,eyeT,varargin)

% Set defaults
print = 0;
nPrint = 500;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'-p','print'},
            print = varargin{varStrInd(iv)+1};
    end
end

inEvTms(inEvs <= 0) = [];
inEvs(inEvs <= 0) = [];
trStartCode = 1666;
trEndCode = 1667;
antiHeader = 1509;
infoEndCode = 2999;

% Find trial starts and ends
startInds = find(inEvs == trStartCode);
endInds = find(inEvs == trEndCode);
antiHeadInds = find(inEvs == antiHeader);
antiHeadInds = antiHeadInds(2:end);
endInfos = find(inEvs == infoEndCode);

if length(antiHeadInds) > length(endInfos),
    while length(antiHeadInds) > length(endInfos),
        antiHeadInds = antiHeadInds(2:end);
    end
elseif length(antiHeadInds) < length(endInfos),
    while length(antiHeadInds) < length(endInfos),
        endInfos = endInfos(2:end);
    end
end

% while (length(antiHeadInds) > length(endInfos)),
%     check=[antiHeadInds((1:length(endInfos))+1),endInfos];
%     antiHeadInds(find(check(:,1) < check(:,2),1)) = [];
% end
    
% Set some constants
stimHs = [.5, 1, 2];
stimVs = [2, 1, .5];

%% Start trial loop
nTrs = length(antiHeadInds)-1;%min([length(startInds),length(endInds)]);
for it = 1:nTrs,
    
    trialCodes{it,1} = inEvs(antiHeadInds(it):endInfos(it));
    trialTimes{it,1} = inEvTms(antiHeadInds(it):endInfos(it));
end
if nTrs < 100,
    Task.Good = 0;
    return;
end
Task.Good = 1;
evntCnt = cellfun(@length,trialCodes);
allCodes = nan(nTrs,max(evntCnt));
% Task = struct('StimOnsetToTrial',[],'SRT',
Task.StimOnsetToTrial = nan(nTrs,1); 
Task.SRT = nan(nTrs,1);
Task.SaccEnd = nan(nTrs,1);
Task.Reward = nan(nTrs,1);
Task.Tone = nan(nTrs,1);
Task.RewardTone = nan(nTrs,1);
Task.ErrorTone = nan(nTrs,1);
Task.FixSpotOn = nan(nTrs,1);
Task.FixSpotOff = nan(nTrs,1);
Task.SingletonColor = nan(nTrs,1);
Task.DistractorColor = nan(nTrs,1);
Task.CueType = nan(nTrs,1);
Task.Eccentricity = nan(nTrs,1);
Task.EndAngle = nan(nTrs,1);
Task.EndEcc = nan(nTrs,1);
Task.SetSize = nan(nTrs,1);
Task.StimLoc = nan(nTrs,4);
Task.StimDiff = nan(nTrs,4);
Task.TargetLoc = nan(nTrs,1);
Task.TrialType = cell(nTrs,1);
Task.OppType = cell(nTrs,1);
Task.Congruent = nan(nTrs,1);
Task.Correct = nan(nTrs,1);
Task.Abort = nan(nTrs,1);
Task.Error = nan(nTrs,1);

for it = 1:nTrs,
    if print && mod(it,nPrint) == 0,
        if exist('printStr','var'),
            for ib = 1:length(printStr),
                fprintf('\b');
            end
        end
        printStr = sprintf('Working on Trial %d of %d...\n',it,nTrs);
        fprintf(printStr);
    end
%     trCodes = inEvs(startInds(it):endInds(it))
    trCodes = inEvs(antiHeadInds(it):endInfos(it));
    trTimes = inEvTms(antiHeadInds(it):endInfos(it));
    allCodes(it,1:evntCnt(it)) = trCodes;
%     try
    %     if any(trCodes == events.Correct_),
    %         Task.
        %% Get timings
        Task.StimOnsetToTrial(it)      = getEvTime(trCodes,trTimes,events.Target_);
    %     Task.refresh_offset(it) = Task.StimOnsetToTrial(it) - tempoStimOn;
        Task.SRT(it)            = getEvTime(trCodes,trTimes,events.Saccade_);
        Task.SaccEnd(it)        = getEvTime(trCodes,trTimes,events.Decide_);
        Task.Reward(it)         = getEvTime(trCodes,trTimes,events.Reward_);
        Task.Tone(it)           = getEvTime(trCodes,trTimes,events.Tone_);    
        Task.RewardTone(it)     = getEvTime(trCodes,trTimes,events.Reward_tone);
        Task.ErrorTone(it)      = getEvTime(trCodes,trTimes,events.Error_tone);
        Task.FixSpotOn(it)      = getEvTime(trCodes,trTimes,events.FixSpotOn_);
        Task.FixSpotOff(it)     = getEvTime(trCodes,trTimes,events.FixSpotOff_);

        %% Get Stim Colors
        colorInfo = unique(trCodes(trCodes >= 700 & trCodes < 800));
        if any(colorInfo >= 700 & colorInfo < 710),
            Task.SingletonColor(it) = colorInfo(colorInfo >= 700 & colorInfo < 710) - 700;
        else
            Task.SingletonColor(it) = nan;
        end
        if any(colorInfo >= 710 & colorInfo < 720),
            Task.DistractorColor(it) = colorInfo(colorInfo >= 710 & colorInfo < 720) - 710;
        else
            Task.DistractorColor(it) = nan;
        end
        if any(colorInfo >= 720 & colorInfo < 730),
            Task.CueType(it) = colorInfo(find(colorInfo >= 720 & colorInfo < 730,1)) - 720;
        else
            Task.CueType(it) = nan;
        end
        if any(colorInfo >= 730 & colorInfo < 740),
            Task.CueColor(it) = colorInfo(find(colorInfo >= 730 & colorInfo < 740,1)) - 730;
        else
            Task.CueColor(it) = nan;
        end

        %% Get Stim Locations
        if any (trCodes >= 900 & trCodes < 1000),
            thisCodes = unique(trCodes(trCodes >= 900 & trCodes < 920));
            Task.Eccentricity(it) = thisCodes(1) - 900;
        else
            Task.Eccentricity(it) = nan;
        end
        locCodes = trCodes(trCodes >= 5000 & trCodes < 6000);
        diffCodes = trCodes(trCodes >= 6000 & trCodes <= 7000);
        if isempty(diffCodes),
            continue
        end
        Task.SetSize(it) = length(diffCodes);
%         isBad = 0;
        nLoops = 0;
        while (nLoops < length(diffCodes)) && length(unique(mod(locCodes,360/Task.SetSize(it)))) > 1,
            uCodes = unique(mod(locCodes,360/Task.SetSize(it)));
            [n,c] = hist(mod(locCodes,360/Task.SetSize(it)),uCodes);
            locCodes(mod(locCodes,360/Task.SetSize(it)) == c(n==min(n))) = [];
            nLoops = nLoops+1;
        end
        if nLoops == length(diffCodes),
            continue;
        end
    %     Task.SetSize(it) = length(locCodes);
        if length(locCodes) > 1,
            for is = 1:Task.SetSize(it),
                Task.StimLoc(it,is) = locCodes(is) - 5000;
                Task.StimDiff(it,is) = trCodes(trCodes >= (6000 + (100*(is))) & trCodes < (6000 + (100*(is)+20))) - (6000 + (100*(is)));
            end
            if any(trCodes >= 800 & trCodes < (800+Task.SetSize(it))),
                Task.TargetLoc(it) = Task.StimLoc(it,trCodes(trCodes >= 800 & trCodes < (800 + Task.SetSize(it)))-799);
            end
        end
        if any(isnan(Task.StimLoc(it,:))) || any(hist(Task.StimLoc(it,:),sort(unique(Task.StimLoc(it,:)))) > 1)
            continue
        end
        if ~isnan(Task.TargetLoc(it)),
            stimInd = Task.StimDiff(it,Task.StimLoc(it,:)==Task.TargetLoc(it))+1;
            if stimInd <= length(stimHs),
                tmpH = stimHs(stimInd);
                tmpV = stimVs(stimInd);
            else
                continue;
            end
            if abs(tmpH - tmpV) < .001,
                Task.TrialType{it} = 'catch';
            elseif tmpH > tmpV,
                Task.TrialType{it} = 'anti';
            elseif tmpV > tmpH,
                Task.TrialType{it} = 'pro';
            else
                keyboard
            end
            clear tmpH tmpV
            
            stimInd = Task.StimDiff(it,Task.StimLoc(it,:)==mod(Task.TargetLoc(it)+180,360))+1;
            if stimInd <= length(stimHs),
                tmpH = stimHs(stimInd);
                tmpV = stimVs(stimInd);
            else
                continue;
            end
            if abs(tmpH - tmpV) < .001,
                Task.OppType{it} = 'catch';
            elseif tmpH > tmpV,
                Task.OppType{it} = 'anti';
            elseif tmpV > tmpH,
                Task.OppType{it} = 'pro';
            else
                keyboard
            end
            clear tmpH tmpV

            if (strcmpi(Task.TrialType{it},'pro') && strcmpi(Task.OppType{it},'anti')) || (strcmpi(Task.TrialType{it},'anti') && strcmpi(Task.OppType{it},'pro')),
                Task.Congruent(it) = 1;
            elseif strcmpi(Task.TrialType{it},Task.OppType{it}),
                Task.Congruent(it) = 0;
            else
                Task.Congruent(it) = 2;
            end
        end

        % Get Trial Outcome
        if any(trCodes == events.Correct_),
            Task.Correct(it) = 1;
        else
            Task.Correct(it) = 0;
        end
        if isnan(Task.SaccEnd(it)) && ~isnan(Task.SRT(it)),
            tmpT = eyeT(eyeT >= Task.StimOnsetToTrial(it) & eyeT <= Task.SRT(it)+3000);
            tmpX = eyeX(eyeT >= Task.StimOnsetToTrial(it) & eyeT <= Task.SRT(it)+3000);
            tmpY = eyeY(eyeT >= Task.StimOnsetToTrial(it) & eyeT <= Task.SRT(it)+3000);
            if strcmpi(Task.TrialType{it},'pro'),
                inWind = klInBox(tmpX,tmpY,[Task.Eccentricity(it)*cos(klDeg2Rad(Task.TargetLoc(it))),Task.Eccentricity(it)*sin(klDeg2Rad(Task.TargetLoc(it)))],4);
            else
                inWind = klInBox(tmpX,tmpY,[Task.Eccentricity(it)*cos(klDeg2Rad(mod(Task.TargetLoc(it)+180,360))),Task.Eccentricity(it)*sin(klDeg2Rad(mod(Task.TargetLoc(it)+180,360)))],4);
            end
            if any(inWind),
                Task.SaccEnd(it) = tmpT(find(inWind,1));
            end
            if ~isnan(Task.SaccEnd(it)),
                if Task.SaccEnd(it) > (100 + Task.SRT(it)),
%                     if Task.Correct(it) == 1 && strcmpi(Task.TrialType{it},'pro'),
%                         keyboard
%                     end
                    Task.Correct(it) = 0;
                end
            end
        end
        
        if any(trCodes == events.Abort_),
            Task.Abort(it) = 1;
        else
            Task.Abort(it) = 0;
        end
        if ~Task.Abort(it) && ~Task.Correct(it),
            Task.Error(it) = 1;
        else
            Task.Error(it) = 0;
        end
        if Task.Abort(it) == 0 && ~isnan(Task.SRT(it)),
            mnX = nanmean(eyeX(eyeT >= (Task.SRT(it) + 50) & eyeT <= (Task.SRT(it) + 100)));
            mnY = nanmean(eyeY(eyeT >= (Task.SRT(it) + 50) & eyeT <= (Task.SRT(it) + 100)));
            Task.EndEcc(it)   = sqrt(mnY^2 + mnX^2);
            baseAngle = klRad2Deg(atan(abs(mnY)/abs(mnX)));
            if (mnX > 0) && (mnY > 0),
                Task.EndAngle(it) = baseAngle;
            elseif (mnX < 0) && (mnY > 0),
                Task.EndAngle(it) = 180-baseAngle;
            elseif (mnX < 0) && (mnY < 0),
                Task.EndAngle(it) = baseAngle + 180;
            elseif (mnX > 0) && (mnY < 0),
                Task.EndAngle(it) = 360-baseAngle;
            end
%             Task.EndAngle(it) = klRad2Deg(asin(mnY/Task.EndEcc(it)));
%             if mnX < 0,
%                 Task.EndAngle(it) = mod(Task.EndAngle(it)+180,360);
%             end
        end
            
%     catch
%     end
    
%     % Get Trial Type
%     if any(trCodes >= 500 & trCodes < 510),
%         isCatch = trCodes(trCodes >= 500 & trCodes < 510)-500;
%     else
%         isCatch = 0;
%     end
%     if any(trCodes >= 600 & trCodes < 610),
%         isAnti = trCodes(trCodes >= 600 & trCodes < 610)-600;
%     else
%         isAnti = 0;
%     end
%     
%     if isCatch,
%         Task.TrialType{it} = 'catch';
%     elseif isAnti,
%         Task.TrialType{it} = 'anti';
%     elseif ~isAnti,
%         Task.TrialType{it} = 'pro';
%     end
    
    
    
    
%     pause
end

%% Get Trial Outcomes
% Correct Trials:
% Task.Correct  = logical(sum(bsxfun(@eq, allCodes, events.Correct_),2)) | logical(sum(bsxfun(@eq,allCodes,events.CatchCorrect_),2));
% 
% % Errors:
% % This needs to be checked and coded more dynamically
% Task.error_names = {'False', 'Early', 'Late', 'FixBreak', 'HoldError', ...
%     'CatchErrorGo', 'CatchErrorNoGo'};
% 
% Task.error = nan(nTrs,1);
% Task.error(Task.Correct == 1) = 0;
% 
% false_resp = logical(sum(bsxfun(@eq, allCodes, events.Error_sacc),2));
% if(any(false_resp))
%     Task.error(false_resp) = 1;
% end
% 
% early_resp = logical(sum(bsxfun(@eq, allCodes, events.EarlySaccade_),2));
% if(any(early_resp))
%     Task.error(early_resp) = 2;
% end
% 
% fix_break = logical(sum(bsxfun(@eq, allCodes, events.FixError_),2));
% if(any(fix_break))
%     Task.error(fix_break) = 4;
% end
% 
% hold_err = logical(sum(bsxfun(@eq, allCodes, events.BreakTFix_),2));
% if(any(hold_err))
%     Task.error(hold_err) = 5;
% end
% 
% catch_go = logical(sum(bsxfun(@eq, allCodes, events.CatchIncorrectG_),2));
% if(any(catch_go))
%     Task.error(catch_go) = 6;
% end
% 
% catch_hold = logical(sum(bsxfun(@eq, allCodes, events.CatchIncorrectG_),2));
% if(any(catch_hold))
%     Task.error(catch_hold) = 7;
% end


Task.AlignTimes = Task.StimOnsetToTrial;
Task.SaccEnd    = Task.SaccEnd          - Task.StimOnsetToTrial;
Task.Reward     = Task.Reward           - Task.StimOnsetToTrial;
Task.Tone       = Task.Tone             - Task.StimOnsetToTrial;
Task.RewardTone = Task.RewardTone       - Task.StimOnsetToTrial;
Task.ErrorTone  = Task.ErrorTone        - Task.StimOnsetToTrial;
Task.FixSpotOn  = Task.FixSpotOn        - Task.StimOnsetToTrial;
Task.FixSpotOff = Task.FixSpotOff       - Task.StimOnsetToTrial;
Task.StimOnset  = Task.StimOnsetToTrial - Task.StimOnsetToTrial; % should be all zero aferwards
Task.GoCue 		= Task.FixSpotOff;
Task.SRT        = Task.SRT              - Task.StimOnsetToTrial - Task.GoCue;

    function outTime = getEvTime(codes,times,match),
        tmp = find(codes==match,1);
        if ~isempty(tmp),
            outTime = times(tmp);
        else
            outTime = nan;
        end
    end
end
    