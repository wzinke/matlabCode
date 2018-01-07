function [allSDF,allSDFTimes,allAreas,allRows,allMonks,allOtherTasks,allSess,allChans,allParams,allTST] = klPullBradenSEF(task,varargin)

switch task
    case {'MG','mg'}
        xlFile = 'klBradenSEFBookKeeping_mg.xlsx';
        sessHead = 'MG';
        headRow = 4;
    case {'search','Search','Capture','capture','cap'}
        xlFile = 'klBradenSEFBookKeeping_search.xlsx';
        sessHead = 'SEARCH';
        headRow = 4;
end

% Set criteria
minRate = 5;
isiCrit = inf;
% snrCrit = .85;
areaCrit = {'FEF','SEF','SC'};
clip = 0;

varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-a','area'}
            areaCrit = varargin{varStrInd(iv)+1};
        case {'-r','reward'}
            doReward = varargin{varStrInd(iv)+1};
        case {'-c'}
            clip = varargin{varStrInd(iv)+1};
    end
end

% Set some constants
monks = {'Fechner','Zorro'};
tebaBase = [tebaMount,'data';

% Set defaults
blWind = -300:-200;
vCheck = 50:150;
mCheck = -50:0;
rCheck = -50:50;
vWind = -300:500;
mWind = -500:300;
rWind = -200:200;

if doReward
    nAligns = 3;
else
    nAligns = 2;
end

allSDF = {};
allSDFTimes = {};
allAreas = {};
allRows = [];
allMonks = {};
allTST = [];
allOtherTasks = {};
allSess = {};
allChans = {};
% Start monk loop
for im = 1:length(monks)
    fprintf('Pulling Braden data from monkey %s...\n',monks{im});
    
    [num,text,all] = xlsread(xlFile,monks{im});
    
    % Get mean rate, ISI, and area columns
    mnRateCol = find(cellfun(@(x) strcmpi(x,'mnRate'),all(headRow,:)),1);
    isiCol = find(cellfun(@(x) strcmpi(x,' < 3ms'),all(headRow,:)),1);
    areaCol = find(cellfun(@(x) strcmpi(x,'Area'),all(headRow,:)),1);
    sessCol = find(cellfun(@(x) strcmpi(x,'Session'),all(headRow,:)),1);
    channelCol = find(cellfun(@(x) strcmpi(x,'Channel'),all(headRow,:)),1);
    unitCol = find(cellfun(@(x) strcmpi(x,'Unit'),all(headRow,:)),1);
%     snrCol = find(cellfun(@(x) strcmpi(x,'SNR'),all(headRow,:)),1);
    
    % Use criteria to pick out good rows
    mnRates = cellfun(@(x) single(x),all((headRow+1):end,mnRateCol),'UniformOutput',0); [mnRates{cellfun(@isempty,mnRates)}] = deal(nan);
    isis = cellfun(@(x) single(x),all((headRow+1):end,isiCol),'UniformOutput',0); [isis{cellfun(@isempty,isis)}] = deal(nan);
%     snrs = cellfun(@(x) single(x),all((headRow+1):end,snrCol),'UniformOutput',0); [snrs{cellfun(@isempty,snrs)}] = deal(nan);
    goodRows = cellfun(@(x) x > minRate,mnRates) & cellfun(@(x) x.*100 < isiCrit,isis) & ismember(all((headRow+1):end,areaCol),areaCrit);% & cellfun(@(x) x >= snrCrit,snrs);
    
    rowInds = find(goodRows)+headRow;
    
    monkSDF = cell(length(rowInds),nAligns);
    monkSDFTimes = cell(length(rowInds),nAligns);
    monkAreas = all(rowInds,areaCol);
    monkMonks = cell(length(rowInds),1); [monkMonks{1:length(rowInds)}] = deal(monks{im});
    monkTST = nan(length(rowInds),1);
    monkOtherTasks = cell(length(rowInds),1);
    monkSess = all(rowInds,sessCol);
    monkChans = cellfun(@(x,y) [num2str(x),num2abc(y)],all(rowInds,channelCol),all(rowInds,unitCol),'UniformOutput',0);
    [monkParams(1:length(rowInds))] = deal(struct('theta',[],'sigma',[],'alpha',[],'offset',[],'blRate',[],'Fano',[],'CV',[],'CV2',[],'LV',[],'LVR',[],'vGOF',[],'mGOF',[],'spkWidth',[]));
    
    % Start row loop
%     for ir = length(rowInds):-1:1
    for ir = 1:length(rowInds)
        myRow = rowInds(ir);
        printHead = sprintf('Analyzing row %d (%d of %d):  ',myRow,ir,length(rowInds));
        fprintf(printHead);

        % Find and load the appropriate .mat files for spikes and behavior
        printStr = 'Loading data...';
        fprintf(printStr);
        
        if (ir > 1) && strcmpi(all{rowInds(ir-1),sessCol},all{rowInds(ir),sessCol})
%         if (ir < length(rowInds)) && strcmpi(all{rowInds(ir+1),sessCol},all{rowInds(ir),sessCol})
        else
            load([tebaBase,filesep,monks{im},filesep,'SEF Popout/Matlab/Final sorts',filesep,all{rowInds(ir),sessCol}]);
        end
        allUnitFiles = dir([tebaBase,filesep,monks{im},filesep,'SEF Popout/Matlab/Final sorts',filesep,strrep(all{rowInds(ir),sessCol},'MG','*')]);
        allTasks = cellfun(@(x) x((strfind(x,'_')+1):(strfind(x,'.')-1)),{allUnitFiles.name},'UniformOutput',0);
        monkOtherTasks{ir} = allTasks(~ismember(allTasks,sessHead));
        fprintf(repmat('\b',1,length(printStr)));
        printStr = 'Calculating SDFs...';
        fprintf(printStr);
        
        Target_(Target_(:,2).*45 > 360,2) = nan;
        
        chanStr = num2str(all{myRow,channelCol}); while length(chanStr) < 2, chanStr = ['0',chanStr]; end
        dspStr = ['DSP',chanStr,num2abc(all{myRow,unitCol})];
        eval(['mySpikes=',dspStr,';']);
        mySpikes(mySpikes == 0) = nan;
        
        % Get the SDFs
        badV = 0; badM = 0; badR = 0;
        spikeMat = mySpikes-repmat(Target_(:,1),1,size(mySpikes,2));
        mMat = spikeMat-repmat(SRT(:,1),1,size(spikeMat,2));
%         rMat = spikeMat-repmat(Reward_(:,3),1,size(spikeMat,2));
        rMat = spikeMat-repmat(SRT(:,1)+500,1,size(spikeMat,2));
        if clip
            spikeMat(spikeMat >= repmat(SRT(:,1),1,size(spikeMat,2))) = nan;
        end
        if sum(isfinite(spikeMat(:))) < 5
            vSDF = 0; vTimes = 0; badV = 1;
        else
            [vSDF,vTimes] = klSpkRatev2(spikeMat,'-q',1);
        end
        if sum(isfinite(mMat(:))) < 5
            mSDF = 0; mTimes = 0; badM = 1;
        else
            [mSDF,mTimes] = klSpkRatev2(mMat,'-q',1);
        end
        if sum(isfinite(rMat(:))) < 5
            rSDF = 0; rTimes = 0; badR = 1;
        else
            [rSDF,rTimes] = klSpkRatev2(rMat,'-q',1);
        end
        
        % Get RFs 
        uLocs = nunique(Target_(:,2));
        vResp = nan(1,length(uLocs));
        mResp = nan(1,length(uLocs));
        rResp = nan(1,length(uLocs));
        for il = 1:length(uLocs)
            vResp(il) = nanmean(nanmean(vSDF(Correct_(:,2)==1 & Target_(:,2)==uLocs(il),ismember(vTimes,vCheck)),1));
            mResp(il) = nanmean(nanmean(mSDF(Correct_(:,2)==1 & Target_(:,2)==uLocs(il),ismember(mTimes,mCheck)),1));
            rResp(il) = nanmean(nanmean(rSDF(Correct_(:,2)==1 & Target_(:,2)==uLocs(il),ismember(rTimes,rCheck)),1));
        end
        if ~any(~isnan(vResp)), vrf = nan; else vrf=uLocs(vResp==max(vResp)); end; if isnan(vrf) || isempty(vrf), badV = 1; end
        if ~any(~isnan(mResp)), mrf = nan; else mrf=uLocs(mResp==max(mResp)); end; if isnan(mrf) || isempty(mrf), badM = 1; end
        if ~any(~isnan(rResp)), rrf = nan; else rrf=uLocs(rResp==max(rResp)); end; if isnan(rrf) || isempty(rrf), badR = 1; end
        
        % Get spike params
        clear vFit mFit
        if badV
            [vFit.mu,vFit.sig,vFit.amp,vFit.bl] = deal(nan);
            vFit.gof = struct('sse',nan,'rsquare',nan,'dfe',nan,'adjrsquare',nan,'rmse',nan);
        else
            vFit = fitGauss(Target_(:,2).*45,Correct_(:,2),vSDF,vTimes,vCheck);
        end
        if badM
            [mFit.mu,mFit.sig,mFit.amp,mFit.bl] = deal(nan);
            mFit.gof = struct('sse',nan,'rsquare',nan,'dfe',nan,'adjrsquare',nan,'rmse',nan);
        else
            mFit = fitGauss(Target_(:,2).*45,Correct_(:,2),mSDF,mTimes,mCheck);
        end
        monkParams(ir).theta = [vFit.mu,mFit.mu];
        monkParams(ir).sigma = [vFit.sig,mFit.sig];
        monkParams(ir).alpha = [vFit.amp,mFit.amp];
        monkParams(ir).offset = [vFit.bl,mFit.bl];
        monkParams(ir).vGOF = vFit.gof;
        monkParams(ir).mGOF = mFit.gof;
        monkParams(ir).blRate = nanmean(nanmean(vSDF(:,ismember(vTimes,blWind)),2),1);
        monkParams(ir).CV = klGetCV(spikeMat);
        monkParams(ir).CV2 = klGetCV(spikeMat,'-type','local');
        monkParams(ir).LV = klGetLV(spikeMat);
        monkParams(ir).LVR = klGetLV(spikeMat,'-type','revised');
        monkParams(ir).Fano = klGetFano(spikeMat);
        monkParams(ir).spkWidth = nan;%klWvWidthv3();
        
        switch task
            case {'MG','mg'}
                if ~badM
                    monkSDF{ir,2} = nanmean(mSDF(ismember(Target_(:,2),mrf) & Correct_(:,2)==1,ismember(mTimes,mWind)),1);
                    monkSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
                else
                    monkSDF{ir,2} = nan(1,length(mWind));
                    monkSDFTimes{ir,2} = mWind;
                end
%                 if ~badV
%                     monkSDF{ir,3} = nanmean(vSDF(ismember(Target_(:,2),vrf) & Correct_(:,2)==1,ismember(vTimes,vWind)),1);
%                     monkSDFTimes{ir,3} = vTimes(ismember(vTimes,vWind));
%                 else
%                     monkSDF{ir,3} = nan(1,length(vWind));
%                     monkSDFTimes{ir,3} = vTimes(ismember(vTimes,vWind));
%                 end
                if doReward
                    if ~badR
                        monkSDF{ir,3} = nanmean(rSDF(ismember(Target_(:,2),rrf) & Correct_(:,2)==1,ismember(rTimes,rWind)),1);
                        if sum(isnan(monkSDF{ir,3}))==length(monkSDF{ir,3})
                            keyboard
                        end
                        monkSDFTimes{ir,3} = rTimes(ismember(rTimes,rWind));
                    else
                        monkSDF{ir,3} = nan(1,length(rWind));
                        monkSDFTimes{ir,3} = rWind;
                    end
                    if isempty(monkSDFTimes{ir,3})
                        keyboard
                    end
                end
                    
                % Clip post-mov spikes from vis to get monkSDF{1}
                clipSpikes = spikeMat;
                clipSpikes(spikeMat > repmat(SRT(:,1),1,size(spikeMat,2))) = nan;
                if sum(isfinite(clipSpikes(:))) < 5
                    monkSDF{ir,1} = nan(1,length(vWind));
                    monkSDFTimes{ir,1} = vWind;
                else
                    [cSDF,cTimes] = klSpkRatev2(clipSpikes,'-q',1);
                    monkSDF{ir,1} = nanmean(cSDF(ismember(Target_(:,2),vrf) & Correct_(:,2)==1,ismember(cTimes,vWind)),1);
                    monkSDFTimes{ir,1} = cTimes(ismember(cTimes,vWind));
                end
            case {'Search','Capture','search','capture','Cap','cap'}
                % Row 1: Target in RF
                % Row 2: Non-Salient Distractor in RF
                % Row 3: Salient Distractor in RF
                % Row 4: Any Distractor in RF
                monkSDF{ir,2} = nan(2,length(mWind));
                monkSDF{ir,1} = nan(2,length(vWind));
                if ~badM
                    monkSDF{ir,2}(1,:) = nanmean(mSDF(ismember(Target_(:,2),mrf) & Correct_(:,2)==1,ismember(mTimes,mWind)),1);
                    %monkSDF{ir,2}(2,:) = nan(1,length(mWind));
                    %monkSDF{ir,2}(3,:) = nanmean(mSDF(~ismember(Target_(:,2),mrf) & Correct_(:,2)==1,ismember(mTimes,mWind)),1);
                    monkSDF{ir,2}(2,:) = nanmean(mSDF(~ismember(Target_(:,2),mrf) & Correct_(:,2)==1,ismember(mTimes,mWind)),1);
                    monkSDFTimes{ir,2} = mTimes(ismember(mTimes,mWind));
                else
                    monkSDF{ir,2} = nan(4,length(mWind));
                    monkSDFTimes{ir,2} = mWind;
                end
                if ~badV
                    % Do running Kruskal-Wallis
                    nWinds = 1;
                    clear p centerWind
                    myWind = -200:-180;
                    centerWind(nWinds) = nanmean(myWind);
                    
                    while sum(ismember(myWind,vTimes))==length(myWind)
                        in1 = nanmean(vSDF(ismember(Target_(:,2),vrf) & Correct_(:,2)==1,ismember(vTimes,myWind)),2);
                        in2 = nanmean(vSDF(~ismember(Target_(:,2),vrf) & Correct_(:,2)==1,ismember(vTimes,myWind)),2);
%                         myNans = isnan(in1) | isnan(in2)
                        if sum(isnan(in1))==length(in1) || sum(isnan(in2))==length(in2)
                            p(nWinds) = nan;
                        else
                            [p(nWinds)] = ranksum(in1,in2);
                        end
                        myWind = myWind+10;
                        nWinds = nWinds+1;
                        centerWind(nWinds)=nanmean(myWind);
                    end
                    nWinds = nWinds-1;
                    centerWind = centerWind(1:(end-1));
                    if any(klGetConsecutive(p < .05) >= 5)
                        tmp = find(klGetConsecutive(p < .05) >= 5 & centerWind >= 50 & centerWind <= (nanmean(SRT(:,1))-50),1);
                        if ~isempty(tmp)
                            monkTST(ir) = tmp;
                        else
                            monkTST(ir) = nan;
                        end
                    else
                        monkTST(ir) = nan;
                    end
                        
                    monkSDF{ir,1} = nan(2,length(vWind));
                    monkSDF{ir,1}(1,1:sum(ismember(vTimes,vWind))) = nanmean(vSDF(ismember(Target_(:,2),vrf) & Correct_(:,2)==1,ismember(vTimes,vWind)),1);
                    %monkSDF{ir,1}(2,:) = nan(1,length(vWind));
                    %monkSDF{ir,1}(3,:) = nanmean(vSDF(~ismember(Target_(:,2),vrf) & Correct_(:,2)==1,ismember(vTimes,vWind)),1);
                    monkSDF{ir,1}(2,1:sum(ismember(vTimes,vWind))) = nanmean(vSDF(~ismember(Target_(:,2),vrf) & Correct_(:,2)==1,ismember(vTimes,vWind)),1);
                    monkSDFTimes{ir,1} = vWind;%vTimes(ismember(vTimes,vWind));
                else
                    monkSDF{ir,1} = nan(4,length(vWind));
                    monkSDFTimes{ir,1} = vWind;
                end
                if doReward
                    if ~badR
%                         monkSDF{ir,3}(1,:) = nanmean(rSDF(ismember(Target_(:,2),rrf) & Correct_(:,2)==1,ismember(rTimes,rWind)),1);
%                         monkSDF{ir,3}(2,:) = nan(1,length(rWind));
%                         monkSDF{ir,3}(3,:) = nanmean(rSDF(~ismember(Target_(:,2),rrf) & Correct_(:,2)==1,ismember(rTimes,rWind)),1);
%                         monkSDF{ir,3}(4,:) = nanmean(rSDF(~ismember(Target_(:,2),rrf) & Correct_(:,2)==1,ismember(rTimes,rWind)),1);
                        monkSDF{ir,3}(1,:) = nanmean(rSDF(Correct_(:,2)==1,ismember(rTimes,rWind)),1);
                        monkSDF{ir,3}(2,:) = nanmean(rSDF(Errors_(:,5)==1,ismember(rTimes,rWind)),1);
                        monkSDFTimes{ir,3} = rTimes(ismember(rTimes,rWind));
                    else
                        monkSDF{ir,3} = nan(4,length(rWind));
                        monkSDFTimes{ir,3} = rWind;
                    end
                end
        end
       fprintf(repmat('\b',1,length(printStr)+length(printHead)));
    end
    allSDF = cat(1,allSDF,monkSDF);
    allSDFTimes = cat(1,allSDFTimes,monkSDFTimes);
    allAreas = cat(1,allAreas,monkAreas);
    allMonks = cat(1,allMonks,monkMonks);
    allRows = cat(1,allRows,rowInds);
    allTST = cat(1,allTST,monkTST);
    allOtherTasks = cat(1,allOtherTasks,monkOtherTasks);
    allSess = cat(1,allSess,monkSess);
    allChans = cat(1,allChans,monkChans);
    if ~exist('allParams','var')
        allParams = monkParams;
    else
        allParams((1:length(monkParams))+length(allParams)) = monkParams;
    end
end 

end


function gFit = fitGauss(Target_,Correct_,sdf,times,wind)

uLocs = nunique(Target_);
% Cut locations that don't make sense.
[n,c] = hist(mod(uLocs,45),nunique(mod(uLocs,45)));
uLocs(mod(uLocs,45)~=c(n==max(n))) = [];
resp = nan(size(uLocs));
for il = 1:length(uLocs)
    resp(il) = nanmean(nanmean(sdf(Target_==uLocs(il) & Correct_==1,ismember(times,wind)),2),1);
end
% Get maximum value
maxLocInd = find(resp==max(resp),1);
shiftAmnt = uLocs(maxLocInd)-180;
blInds = mod((-1:1)+mod(maxLocInd+length(resp)/2,length(resp)),length(resp)); blInds(blInds==0) = length(resp);
blSub = nanmean(resp(blInds));
fitX = mod(uLocs(~isnan(resp)),360);
fitY = resp(~isnan(resp))-blSub;
if sum(isfinite(fitY)) < 3
    [gFit.amp,gFit.mu,gFit.sig,gFit.bl] = deal(nan);
    gFit.gof = struct('sse',nan,'rsquare',nan,'dfe',nan,'adjrsquare',nan,'rmse',nan);
    return
end
try
    [f, gof]=fit(fitX,fitY,'gauss1');
catch
    [f.a1, f.b1, f.c1] = deal(nan);
    gof = struct('sse',nan,'rsquare',nan,'dfe',nan,'adjrsquare',nan,'rmse',nan);
end
gFit.amp = f.a1;
gFit.mu = f.b1+shiftAmnt;
gFit.sig = f.c1;
gFit.bl = blSub;
gFit.gof = gof;
end
        
        
        
        