function Task = klAnalyzeAnti(fileName,varargin),

% Set defaults
print = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)}
        case {'-p','print'},
            print = varargin{varStrInd(iv)+1};
    end
end

if print,
    printStr = 'Getting Events...';
    fprintf(printStr);
end

tdtEvRaw = TDT2mat(fileName,'TYPE',{'epocs'},'VERBOSE',0);
if isfield(tdtEvRaw.epocs,'STRB'),
    tdtEvs = tdtEvRaw.epocs.STRB.data;
    tdtEvTms = tdtEvRaw.epocs.STRB.onset.*1000;
else
    tdtEvRaw = TDT2mat(fileName,'TYPE',{'scalars'},'VERBOSE',0);
    tdtEvs = tdtEvRaw.scalars.EVNT.data;
    if any(tdtEvs >= (2^15)),
        tdtEvs = tdt2EvShft(tdtEvs);
    end
    if any(mod(tdtEvs,1)) ~= 0,
        tdtEvs = tdtEvRaw.scalars.EVNT.data - (2^15);
    end
    tdtEvTms = tdtEvRaw.scalars.EVNT.ts.*1000;
    tdtEvTms(tdtEvs < 0) = [];
    tdtEvs(tdtEvs < 0) = [];
end
EV = TEMPO_EV_cosman_rig028;
if print,
    for ib = 1:length(printStr),
        fprintf('\b');
    end
    printStr = 'Getting Eye X...';
    fprintf(printStr);
end
eyeXRaw = TDT2mat(fileName,'TYPE',{'streams'},'STORE','EyeX','VERBOSE',0);
if print,
    for ib = 1:length(printStr),
        fprintf('\b');
    end
    printStr = 'Getting Eye Y...';
    fprintf(printStr);
end
eyeYRaw = TDT2mat(fileName,'TYPE',{'streams'},'STORE','EyeY','VERBOSE',0);
eyeX = eyeYRaw.streams.EyeY.data.*(3);
eyeY = eyeXRaw.streams.EyeX.data.*(-3);
eyeT = (0:(length(eyeX)-1)).*(1000/eyeXRaw.streams.EyeX.fs);
eyeR = sqrt(eyeX.^2 + eyeY.^2);
eyeThRaw = klRad2Deg(atan(abs(eyeY)./abs(eyeX)));
eyeTh = nan(size(eyeThRaw));
eyeTh(eyeX > 0 & eyeY > 0) = eyeThRaw(eyeX > 0 & eyeY > 0);
eyeTh(eyeX < 0 & eyeY > 0) = 180-eyeThRaw(eyeX < 0 & eyeY > 0);
eyeTh(eyeX < 0 & eyeY < 0) = 180+eyeThRaw(eyeX < 0 & eyeY < 0);
eyeTh(eyeX > 0 & eyeY < 0) = 360-eyeThRaw(eyeX > 0 & eyeY < 0);
Eyes.X = eyeX;
Eyes.Y = eyeY;
Eyes.R = eyeR;
Eyes.Theta = eyeTh;
Eyes.Times = eyeT;
if print,
    for ib = 1:length(printStr),
        fprintf('\b');
    end
    printStr = 'Saving Eyes...';
    fprintf(printStr);
end
save(sprintf('%s/Eyes.mat',fileName),'Eyes','-v7.3');
clear Eyes eyeR eyeThRaw eyeTh
if print,
    for ib = 1:length(printStr),
        fprintf('\b');
    end
    printStr = 'Decoding Event Codes...\n';
    fprintf(printStr);
end
Task = klGetAntiTaskv2(tdtEvs,tdtEvTms,EV,eyeX,eyeY,eyeT,'-p',print);
if print,
%     for ib = 1:length(printStr),
%         fprintf('\b');
%     end
    printStr = 'Saving Behavior...';
    fprintf(printStr);
end
save(sprintf('%s/Behav.mat',fileName),'Task');


if ~Task.Good
    fprintf('\n');
    return;
end

% Some constants...
nTrialsRun = 10;

% Get some handy logical vectors
isPro = strcmpi(Task.TrialType,'pro');
isAnti = strcmpi(Task.TrialType,'anti');
isCong = Task.Congruent == 1;
isIncong = Task.Congruent == 0;
congCatch = Task.Congruent == 2;
isGood = Task.Abort == 0;

% Get all RT vector
allRTs = Task.SRT;% + Task.GoCue;

% Get overall performance for Anti and anti
percPro = sum(Task.Correct == 1 & isGood & isPro)/sum(isGood & isPro);
percAnti = sum(Task.Correct == 1 & isGood & isAnti)/sum(isGood & isAnti);
rtPro = nanmean(allRTs(Task.Correct == 1 & isGood & isPro));
rtAnti = nanmean(allRTs(Task.Correct == 1 & isGood & isAnti));


% Get Congruent/incongruent performance
percProCong = sum(Task.Correct == 1 & isGood & isPro & isCong)/sum(isGood & isPro & isCong);
percProIncong = sum(Task.Correct == 1 & isGood & isPro & isIncong)/sum(isGood & isPro & isIncong);
percProCatch = sum(Task.Correct == 1 & isGood & isPro & congCatch)/sum(isGood & isPro & congCatch);

percAntiCong = sum(Task.Correct == 1 & isGood & isAnti & isCong)/sum(isGood & isAnti & isCong);
percAntiIncong = sum(Task.Correct == 1 & isGood & isAnti & isIncong)/sum(isGood & isAnti & isIncong);
percAntiCatch = sum(Task.Correct == 1 & isGood & isAnti & congCatch)/sum(isGood & isAnti & congCatch);

% Get congruent/incongruent RTs
rtProCong = nanmean(allRTs(Task.Correct == 1 & isGood & isPro & isCong));
rtProIncong = nanmean(allRTs(Task.Correct == 1 & isGood & isPro & isIncong));
rtProCatch = nanmean(allRTs(Task.Correct == 1 & isGood & isPro & congCatch));

rtAntiCong = nanmean(allRTs(Task.Correct == 1 & isGood & isAnti & isCong));
rtAntiIncong = nanmean(allRTs(Task.Correct == 1 & isGood & isAnti & isIncong));
rtAntiCatch = nanmean(allRTs(Task.Correct == 1 & isGood & isAnti & congCatch));

% Now split by location...
uLocs = unique(Task.TargetLoc(~isnan(Task.TargetLoc)));
% rows = locations, columns = [overall, congruent, incongruent, catch]
proPerf = nan(length(uLocs),4);
antiPerf = nan(length(uLocs),4);
locProRTs = nan(length(uLocs),4);
locAntiRTs = nan(length(uLocs),4);

for il = 1:length(uLocs),
    basicLog = isGood & (Task.TargetLoc == uLocs(il));
    
    proPerf(il,1) = sum(isPro & basicLog & Task.Correct == 1)/sum(isPro & basicLog);
    proPerf(il,2) = sum(isPro & basicLog & Task.Correct == 1 & isCong)/sum(isPro & basicLog & isCong);
    proPerf(il,3) = sum(isPro & basicLog & Task.Correct == 1 & isIncong)/sum(isPro & basicLog & isIncong);
    proPerf(il,4) = sum(isPro & basicLog & Task.Correct == 1 & congCatch)/sum(isPro & basicLog & congCatch);
    
    antiPerf(il,1) = sum(isAnti & basicLog & Task.Correct == 1)/sum(isAnti & basicLog);
    antiPerf(il,2) = sum(isAnti & basicLog & Task.Correct == 1 & isCong)/sum(isAnti & basicLog & isCong);
    antiPerf(il,3) = sum(isAnti & basicLog & Task.Correct == 1 & isIncong)/sum(isAnti & basicLog & isIncong);
    antiPerf(il,4) = sum(isAnti & basicLog & Task.Correct == 1 & congCatch)/sum(isAnti & basicLog & congCatch);
    
    locProRTs(il,1) = nanmean(allRTs(isPro & basicLog & Task.Correct == 1));
    locProRTs(il,2) = nanmean(allRTs(isPro & basicLog & Task.Correct == 1 & isCong));
    locProRTs(il,3) = nanmean(allRTs(isPro & basicLog & Task.Correct == 1 & isIncong));
    locProRTs(il,4) = nanmean(allRTs(isPro & basicLog & Task.Correct == 1 & congCatch));

    locAntiRTs(il,1) = nanmean(allRTs(isAnti & basicLog & Task.Correct == 1));
    locAntiRTs(il,2) = nanmean(allRTs(isAnti & basicLog & Task.Correct == 1 & isCong));
    locAntiRTs(il,3) = nanmean(allRTs(isAnti & basicLog & Task.Correct == 1 & isIncong));
    locAntiRTs(il,4) = nanmean(allRTs(isAnti & basicLog & Task.Correct == 1 & congCatch));
end

%% Plot summaries
figure(); bar(uLocs,proPerf(:,2:end))
title('Pro Performance','fontsize',18);
xlabel('Singleton Location','fontsize',16);
ylabel('Percent Correct','fontsize',16);
set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',14,'YTick',0:.25:1,'YTickLabel',0:25:100,'XLim',[min(uLocs)-nanmean(diff(uLocs)),max(uLocs)+2*nanmean(diff(uLocs))]);
l=legend({'Congruent','Incongruent','Square'});
set(gcf,'paperposition',[.2 .1 10.5 7.5]);
saveas(gcf,sprintf('%s/ProPerformance-ByLocation.png',fileName));

figure(); bar(uLocs,antiPerf(:,2:end))
title('Anti Performance','fontsize',18);
xlabel('Singleton Location','fontsize',16);
ylabel('Percent Correct','fontsize',16);
set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',14,'YTick',0:.25:1,'YTickLabel',0:25:100,'XLim',[min(uLocs)-nanmean(diff(uLocs)),max(uLocs)+2*nanmean(diff(uLocs))]);
l=legend({'Congruent','Incongruent','Square'});
set(gcf,'paperposition',[.2 .1 10.5 7.5]);
saveas(gcf,sprintf('%s/AntiPerformance-ByLocation.png',fileName));

figure(); bar(uLocs,locProRTs(:,2:end)); 
title('Pro Reaction Time','fontsize',18);
xlabel('Singleton Location','fontsize',16);
ylabel('Correct RT (ms)','fontsize',16);
set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',14,'XLim',[min(uLocs)-nanmean(diff(uLocs)),max(uLocs)+2*nanmean(diff(uLocs))]);
l=legend({'Congruent','Incongruent','Square'});
set(gcf,'paperposition',[.2 .1 10.5 7.5]);
saveas(gcf,sprintf('%s/ProRTs-ByLocation.png',fileName));

figure(); bar(uLocs,locAntiRTs(:,2:end)); title('RT Anti');
title('Anti Reaction Time','fontsize',18);
xlabel('Singleton Location','fontsize',16);
ylabel('Correct RT (ms)','fontsize',16);
set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',14,'XLim',[min(uLocs)-nanmean(diff(uLocs)),max(uLocs)+2*nanmean(diff(uLocs))]);
l=legend({'Congruent','Incongruent','Square'});
set(gcf,'paperposition',[.2 .1 10.5 7.5]);
saveas(gcf,sprintf('%s/AntiRTs-ByLocation.png',fileName));

%% Now plot time series
% % Start with RTs
% figure();
% plot(find(isPro & isCong & Task.Correct == 1),allRTs(isPro & isCong & Task.Correct == 1),'color','g','linewidth',2,'linestyle','-'); hold on;
% plot(find(isPro & isIncong & Task.Correct == 1),allRTs(isPro & isIncong & Task.Correct == 1),'color','g','linewidth',2,'linestyle','--'); hold on;
% plot(find(isPro & congCatch & Task.Correct == 1),allRTs(isPro & congCatch & Task.Correct == 1),'color','g','linewidth',2,'linestyle',':'); hold on;
% plot(find(isAnti & isCong & Task.Correct == 1),allRTs(isAnti & isCong & Task.Correct == 1),'color','r','linewidth',2,'linestyle','-'); hold on;
% plot(find(isAnti & isIncong & Task.Correct == 1),allRTs(isAnti & isIncong & Task.Correct == 1),'color','r','linewidth',2,'linestyle','--'); hold on;
% plot(find(isAnti & congCatch & Task.Correct == 1),allRTs(isAnti & congCatch & Task.Correct == 1),'color','r','linewidth',2,'linestyle',':'); hold on;

% Now do running average performance
proAnti = [isPro,isAnti];
cong = [isCong,isIncong,congCatch];
% gKern = klMakeGauss(nTrialsRun);
colors = 'gr';
lstyle = {'-','--',':'};

figure();  hold on;
for ic = 1:size(cong,2),
    for ip = 1:size(proAnti,2),
        myTrials = find(isGood & proAnti(:,ip) & cong(:,ic));
        gKern = klMakeGauss(round(length(myTrials)*.05));
        myPerf = Task.Correct(myTrials)';
        convPerf = conv2(myPerf,gKern./sum(gKern),'same');
        plot(myTrials,convPerf.*100,'color',colors(ip),'linestyle',lstyle{ic},'linewidth',2);
    end
    gKern = klMakeGauss(round(sum(isGood & cong(:,2))*.05));
    myTrials = find(isGood & cong(:,ic));
    myPerf = Task.Correct(myTrials)';
    convPerf = conv2(myPerf,gKern./sum(gKern),'same');
    plot(myTrials,convPerf.*100,'color','k','linestyle',lstyle{ic},'linewidth',2);
end
% gKern = klMakeGauss(round(sum(isGood)*.05));
% convPerf = conv2(Task.Correct(isGood),gKern./sum(gKern),'same');
% plot(find(isGood),convPerf.*100,'color','k','linestyle','-','linewidth',2);
set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',14,'YTick',0:25:100);
xlabel('Trial #','fontsize',16);
ylabel('Percent Correct','fontsize',16);
title('Running Average Performance','fontsize',18);
set(gcf,'paperposition',[.2 .1 10.5 7.5]);
saveas(gcf,sprintf('%s/RunningAveragePerformance.png',fileName));

% Running average RTs...
figure();  hold on;
for ic = 1:size(cong,2),
    for ip = 1:size(proAnti,2),
        myTrials = find(isGood & proAnti(:,ip) & cong(:,ic) & Task.Correct == 1);
        if length(myTrials) < 15,
            continue
        end
        myRTs = allRTs(myTrials);
        mm = movmean(myRTs,round(length(myTrials)*.05));%nTrialsRun);
        plot(myTrials,mm,'color',colors(ip),'linestyle',lstyle{ic},'linewidth',2);
    end
    if sum(Task.Correct == 1 & cong(:,ic)) >= 15,
        plot(find(Task.Correct == 1 & cong(:,ic)),movmean(allRTs(Task.Correct == 1 & cong(:,ic)),round(sum(Task.Correct == 1 & cong(:,ic))*.05)),'color','k','linestyle',lstyle{ic},'linewidth',2);
    end
end
% plot(find(Task.Correct' == 1),movmean(allRTs(Task.Correct' == 1),round(sum(Task.Correct == 1)*.05)),'color','k','linestyle','-','linewidth',2);
set(gca,'box','off','tickdir','out','ticklength',get(gca,'ticklength').*3,'fontsize',14);
xlabel('Trial #','fontsize',16);
ylabel('Correct RT (ms)','fontsize',16);
title('Running Average Reaction Times','fontsize',18);
set(gcf,'paperposition',[.2 .1 10.5 7.5]);
saveas(gcf,sprintf('%s/RunningAverageRTs.png',fileName));

%% Let's look at endpoint errors...

% Start with Pro...
[n,c] = hist(Task.SetSize,2:2:12);
if (c(n==max(n))==4),
    % Start with Pro...
    figure();
    uLocs = unique(Task.TargetLoc(~isnan(Task.TargetLoc) & Task.SetSize == 4));
    % In case of a switch, let's grab the top 4 most popular locations
    [n,c] = hist(Task.TargetLoc,uLocs);
    [sortN,sortInds] = sort(n,'descend');
    uLocs = c(sortInds(1:4));
    for il = 1:length(uLocs),
        subplot(2,2,il);
        polarscatter(klDeg2Rad(Task.EndAngle(proAnti(:,1) & ~isnan(Task.EndAngle) & Task.Correct == 1 & Task.TargetLoc==uLocs(il))),Task.EndEcc(proAnti(:,1) & ~isnan(Task.EndAngle) & Task.Correct == 1 & Task.TargetLoc == uLocs(il)),[],'g'); hold all;
        polarscatter(klDeg2Rad(Task.EndAngle(proAnti(:,1) & ~isnan(Task.EndAngle) & Task.Correct == 0 & Task.TargetLoc==uLocs(il))),Task.EndEcc(proAnti(:,1) & ~isnan(Task.EndAngle) & Task.Correct == 0 & Task.TargetLoc == uLocs(il)),[],'r');
        title(sprintf('Singleton Loc = %d',uLocs(il)));
        set(gca,'RLim',[0,max(Task.Eccentricity)+2]);
    end
    set(gcf,'paperposition',[.2 .1 10.5 7.5]);
    saveas(gcf,sprintf('%s/EndPointErrors-Pro-ByLocation.png',fileName));
   
    % Now Anti...
    figure();
    uLocs = unique(Task.TargetLoc(~isnan(Task.TargetLoc) & Task.SetSize == 4));
    % In case of a switch, let's grab the top 4 most popular locations
    [n,c] = hist(Task.TargetLoc,uLocs);
    [sortN,sortInds] = sort(n,'descend');
    uLocs = c(sortInds(1:4));
    for il = 1:length(uLocs),
        subplot(2,2,il);
        polarscatter(klDeg2Rad(Task.EndAngle(proAnti(:,2) & ~isnan(Task.EndAngle) & Task.Correct == 1 & Task.TargetLoc==uLocs(il))),Task.EndEcc(proAnti(:,2) & ~isnan(Task.EndAngle) & Task.Correct == 1 & Task.TargetLoc == uLocs(il)),[],'g'); hold on;
        polarscatter(klDeg2Rad(Task.EndAngle(proAnti(:,2) & ~isnan(Task.EndAngle) & Task.Correct == 0 & Task.TargetLoc==uLocs(il))),Task.EndEcc(proAnti(:,2) & ~isnan(Task.EndAngle) & Task.Correct == 0 & Task.TargetLoc == uLocs(il)),[],'r');
        title(sprintf('Singleton Loc = %d',uLocs(il)));
        set(gca,'RLim',[0,max(Task.Eccentricity)+2]);
    end
    set(gcf,'paperposition',[.2 .1 10.5 7.5]);
    saveas(gcf,sprintf('%s/EndPointErrors-Anti-ByLocation.png',fileName));
   
end

% Regardless of Set Size, let's see about congruency

% First, rotate the end angle by the target location...
rotAngle = Task.EndAngle - Task.TargetLoc;
rotAngle(rotAngle < 0) = rotAngle(rotAngle < 0) + 360;

figure();
colors = 'rcm';
congLabs = {'Cong','Incong','Square'};
for ic = 1:3,
    subplot(3,2,1+2*(ic-1));
%     subplot(1,2,1);
    polarscatter(klDeg2Rad(rotAngle(proAnti(:,1) & Task.Correct == 1 & cong(:,ic))),Task.EndEcc(proAnti(:,1) & Task.Correct == 1 & cong(:,ic)),[],'g'); hold on;
    polarscatter(klDeg2Rad(rotAngle(proAnti(:,1) & Task.Correct == 0 & isGood & cong(:,ic))),Task.EndEcc(proAnti(:,1) & Task.Correct == 0 & isGood & cong(:,ic)),[],colors(ic));
%     polarscatter(klDeg2Rad(rotAngle(proAnti(:,1) & Task.Correct == 0 & isGood & cong(:,2))),Task.EndEcc(proAnti(:,1) & Task.Correct == 0 & isGood & cong(:,2)),[],'c');
%     polarscatter(klDeg2Rad(rotAngle(proAnti(:,1) & Task.Correct == 0 & isGood & cong(:,3))),Task.EndEcc(proAnti(:,1) & Task.Correct == 0 & isGood & cong(:,3)),[],'m');
    if ic == 1,
        title('Pro Trials');
    end
    set(gca,'RLim',[0,max(Task.Eccentricity)+2]);
%     ylabel(congLabs{ic});
    
    subplot(3,2,2+2*(ic-1));
    polarscatter(klDeg2Rad(rotAngle(proAnti(:,2) & Task.Correct == 1 & cong(:,ic))),Task.EndEcc(proAnti(:,2) & Task.Correct == 1 & cong(:,ic)),[],'g'); hold on;
    polarscatter(klDeg2Rad(rotAngle(proAnti(:,2) & Task.Correct == 0 & isGood & cong(:,ic))),Task.EndEcc(proAnti(:,2) & Task.Correct == 0 & isGood & cong(:,ic)),[],colors(ic));
%     polarscatter(klDeg2Rad(rotAngle(proAnti(:,2) & Task.Correct == 0 & isGood & cong(:,2))),Task.EndEcc(proAnti(:,2) & Task.Correct == 0 & isGood & cong(:,2)),[],'c');
%     polarscatter(klDeg2Rad(rotAngle(proAnti(:,2) & Task.Correct == 0 & isGood & cong(:,3))),Task.EndEcc(proAnti(:,2) & Task.Correct == 0 & isGood & cong(:,3)),[],'m');
    if ic == 1,
        title('Anti Trials');
    end
    set(gca,'RLim',[0,max(Task.Eccentricity)+2]);
end
set(gcf,'paperposition',[.2 .1 10.5 7.5]);
saveas(gcf,sprintf('%s/EndPointErrors-Congruency.png',fileName));
   






    