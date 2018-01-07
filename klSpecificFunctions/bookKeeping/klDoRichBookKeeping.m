function klDoRichBookKeeping(task,varargin)

switch task
    case {'MG','mg'}
        outFile = './klRichBookKeeping_mg.xlsx';
        sessHead = 'MG';
    case {'Search','Capture','search','capture','cap','SAT','sat'}
        outFile = './klRichBookKeeping_search.xlsx';
        sessHead = 'SEARCH';
end

monks = {'Darwin','Euler','Quincy','Seymour'};
tebaBase = [tebaMount,'data'];

% Set defaults
blWind = -300:-100;
vWind = -200:300;
mWind = -300:200;
minTST = 0;
alph = .05;
consecCrit = 20;
bsxTarg = 0; bsxDist = -1:1;

% Start monkey loop
for im = 1:length(monks)
    outMat = {};
    
    % Get list of sessions - found in SAT
    satSess = dir([tebaBase,filesep,monks{im},filesep,'SAT/Matlab',filesep,upper(monks{im}(1)),'*',sessHead,'*']);
    satSefSess = dir([tebaBase,filesep,monks{im},filesep,'SAT_SEF/Matlab',filesep,upper(monks{im}(1)),'*',sessHead,'*']);
    allSess = [satSess;satSefSess];
    keepVars = whos; keepNames = {keepVars.name}; keepNames{end+1} = 'keepNames';
    fprintf('Analyzing Rich data for %s...(%d sessions)\n',monks{im},length(allSess));
    
    % Start session loop
    for is = 1:length(allSess)
        fprintf('\tSession %s (%d of %d)...',allSess(is).name,is,length(allSess));
        % Load up the data
        load([allSess(is).folder,filesep,allSess(is).name]);
        
        % Get units
        sessUnits = whos('DSP*');
        switch sessHead
            case {'SEARCH'}
                otherObj = matfile([allSess(is).folder,filesep,allSess(is).name(1:(end-10)),'MG.mat']);
                otherDSPs = whos(otherObj,'DSP*');
                otherDSPs = {otherDSPs.name};
            case {'MG'}
                otherObj = matfile([allSess(is).folder,filesep,allSess(is).name(1:(end-6)),'SEARCH.mat']);
                otherDSPs = whos(otherObj,'DSP*');
                otherDSPs = {otherDSPs.name};
        end
        printStr = [];
        for iu = 1:length(sessUnits)
            fprintf(repmat('\b',1,length(printStr)));
            printStr = sessUnits(iu).name;
            fprintf(printStr);
            
            thisChan = str2double(sessUnits(iu).name(4:5));
            thisUnit = abc2num(sessUnits(iu).name(6:end));
            hasOther = any(ismember(otherDSPs,sessUnits(iu).name));
            eval(['mySpikes=',sessUnits(iu).name,';']);
            mySpikes(mySpikes==0) = nan;
            
            % Get SDFs
            [vSDF,vTimes] = klSpkRatev2(mySpikes-repmat(Target_(:,1),1,size(mySpikes,2)),'-q',1);
            [mSDF,mTimes] = klSpkRatev2(mySpikes-repmat(Target_(:,1)+SRT(:,1),1,size(mySpikes,2)),'-q',1);
            
            % Get baseline rate
            blRate = nanmean(nanmean(vSDF(:,ismember(vTimes,blWind)),2));
            % Get overall mean rate
            catResp = [vSDF(:,ismember(vTimes,vWind)),mSDF(:,ismember(mTimes,mWind))];
            mnRate = nanmean(nanmean(catResp,2),1);
            
            % Get ISI < 3ms
            isiMat = diff(mySpikes,[],2);
            percISI = sum(isiMat(:) < 3)./sum(isfinite(isiMat(:))); 
            
            % Get tuning stuff
            uLocs = nunique(Target_(:,2));
            vResp = nan(length(uLocs),1);
            mResp = nan(length(uLocs),1);
            for il = 1:length(uLocs)
                vResp(il) = nanmean(nanmean(vSDF(Target_(:,2)==uLocs(il) & Correct_(:,2)==1,vTimes >= 50 & vTimes <= nanmean(SRT(:,1))),2),1);
                mResp(il) = nanmean(nanmean(mSDF(Target_(:,2)==uLocs(il) & Correct_(:,2)==1,mTimes >= -100 & mTimes <= 100),2),1);
            end
            if ~any(~isnan(vResp)), vrf = nan; else vrf=uLocs(find(vResp==max(vResp),1)); end; if isnan(vrf) || isempty(vrf), badV = 1; end
            if ~any(~isnan(mResp)), mrf = nan; else mrf=uLocs(find(mResp==max(mResp),1)); end; if isnan(mrf) || isempty(mrf), badM = 1; end

            rayV=circ_rtest(klDeg2Rad(uLocs.*45),vResp);
            rayM = circ_rtest(klDeg2Rad(uLocs.*45),mResp);
%             vRF = RFs.(sessUnits(iu).name);
%             mRF = MFs.(sessUnits(iu).name);
%             
%             [~,whichFold] = fileparts(fileparts(allSess(is).folder));
            avSDF = klRunningAvv2(vSDF,5);
            myTimes = find(vTimes >= -50 & vTimes <= nanmean(SRT(:,1))-50);

            targIn = avSDF(Correct_(:,2)==1 & ismember(Target_(:,2),bsxfun(@plus,vrf,bsxTarg)),myTimes);
            targOut = avSDF(Correct_(:,2)==1 & ~ismember(Target_(:,2),bsxfun(@plus,vrf,bsxDist)),myTimes);
            subTimes = vTimes(myTimes(1:5:end));
            [auc,~,bootP,~] = klROCv2(targIn(:,1:5:length(myTimes)),targOut(:,1:5:length(myTimes)));
            if isnan(auc)
                mnAUC = nan;
                tst = nan;
            else
                mnAUC = nanmean(auc(subTimes >= 50));
                consecP = klGetConsecutive(bootP' < alph);
                if any(consecP >= consecCrit/5)
                    tst = subTimes(find(consecP >= consecCrit/5,1));
                    if tst <= minTST, tst = nan; end
                else
                    tst = nan;
                end
            end
            
%             vRF = RFs.(sessUnits(iu).name);
%             mRF = MFs.(sessUnits(iu).name);
%             
            [~,whichFold] = fileparts(fileparts(allSess(is).folder));
            outRow = {whichFold,allSess(is).name,thisChan,thisUnit,hasOther,blRate,mnRate,percISI,rayV,rayM,vrf,mrf,mnAUC,tst};
            outMat = cat(1,outMat,outRow);
        end
        currVars = whos; currNames = {currVars.name}; keepNames{end+1} = 'currNames';
        currNames = currNames(~ismember(currNames,keepNames));
        for iv = 1:length(currNames)
            eval(['clear ',currNames{iv}]);
        end
        fprintf('\n');
    end
    xlwrite(outFile,outMat,monks{im},'A5');
end 