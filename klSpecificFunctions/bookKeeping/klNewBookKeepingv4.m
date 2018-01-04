%% klBookKeepingNew
clearvars ; close all;
% if ~exist('./klNewBookKeepingv2.m','file'),
%     fprintf('\n*** Please change directory to main MATLAB folder...\n');
%     keyboard
% end

freshRun = 1;

% Initialize some variables
monk = {'Gauss','Helmholtz','Darwin'};
% monk = {'Gauss'};
xlFile = './plexBookKeeping/klDataBookKeeping_mg.xlsx';
dataFold = 'Y:/Users/Wolf/ephys_db';

type        = 'DSP';                             % DSP for spikes, LFP for LFP, and EEG for EEG
task        = 'MG';
physType    = 'rates';
visualize   = 0;
doBL        = 0;
bWind       = 200;
blWind      = [-300,-100];
vWind       = [50 150];
mWind       = [0 50];

% Set filters for nonsense values
rateLim     = 100;
visMax      = 200;
movMin      = -200;
alph        = .1;b

rOff = 4;
cOff = 3;

startCell = [num2abc(cOff+1),num2str(rOff+1)];
makeColLookup(1)
load xlCols

startTic = tic;
% load('sortedWaves.mat');
load meanSorts.mat

for im = 1:length(monk)
    monkTic = tic;
    thisMonk =monk{im};
    
	if ~freshRun,
		[excelNum,~,excelAll] = xlsread(xlFile,monk{im});
	end
	
    % Initialize Variables
    outMat = {};
    startRowOut = 1;
    
    %% Start Loops
    allSess = dir(sprintf('%s/%s/2*',dataFold,thisMonk));
    sessNames = {allSess.name};
    myRow = 0;
    for is = 1:length(sessNames),
        sessTic = tic;
        sessUnits = dir(sprintf('%s/%s/%s/%s/*DSP*',dataFold,thisMonk,sessNames{is},type));
        unitNames = {sessUnits.name};
        if is == 1,
            fprintf('Analyzing session %s (%d/%d - %.1f%%)',sessNames{is},is,length(sessNames),(100*is)/length(sessNames));
        else
            for ib = 1:numBack
                fprintf('\b');
            end
            fprintf('%s (%d/%d - %.1f%%)',sessNames{is},is,length(sessNames),(100*is)/length(sessNames));
        end
        numBack = length(sprintf('%s (%d/%d - %.1f%%)',sessNames{is},is,length(sessNames),(100*is)/length(sessNames)));
        
        chanLoaded = 0;
        for iu = 1:length(unitNames),
            clear f fName
            outRow = cell(1,40);
            
            % Check for valid data files
            f = dir(sprintf('%s/%s/%s/%s/%s/%s_%s_%s*',dataFold,monk{im},sessNames{is},type,unitNames{iu},sessNames{is},unitNames{iu},task));
            if length(f) == 0,
%                 fprintf('\tNo "%s" file found for unit %s... Skipping to next unit\n',task,unitNames{iu});
                continue;
            elseif length(f) > 1,
%                 fprintf('\tMultiple "%s" files found for unit %s... Please check directory structure\n',task,unitNames{iu});
                %keyboard
            end
            fName = f(1).name;
            myRow = myRow+1;
             
			if ~freshRun,
				if any(strcmp(excelAll(:,strcmp(excelAll(4,:),'Name')),sessNames{is}) & strcmp(excelAll(:,strcmp(excelAll(4,:),'chanCode')),unitNames{iu})),
					continue;
                end
                startRowOut = iu; 
			end
			
            % Load unit data
            loadTic = tic;
            
            load(sprintf('%s/%s/%s/%s/%s/%s',dataFold,monk{im},sessNames{is},type,unitNames{iu},fName));
            load(sprintf('%s/%s/%s/%s/%s/waves/%s_waves.mat',dataFold,monk{im},sessNames{is},type,unitNames{iu},fName(1:(end-4))));
            
            if size(spiketimes,2) < 1,
                continue;
            end
            
            [vSDF,vTimes] = klSpkRatev2(spiketimes,'-q',1);
            [mSDF,mTimes] = klSpkRatev2(spiketimes-repmat(Task.SRT+Task.GoCue,1,size(spiketimes,2)),'-q',1);
            blMean = nanmean(nanmean(vSDF(:,vTimes >= blWind(1) & vTimes < blWind(2)),1),2);
            blStd = nanstd(nanmean(vSDF(:,vTimes >= blWind(1) & vTimes < blWind(2)),1),[],2);
            visMean = nanmean(nanmean(vSDF(:,vTimes >= vWind(1) & vTimes <= vWind(2)),1),2);
            movMean = nanmean(nanmean(mSDF(:,mTimes <= -mWind(1) & mTimes >= -mWind(2)),1),2);
            
            %% Get Identification Info
            outRow{col.sess}  =  sessNames{is};
            chanLoaded = chanLoaded + 1;
            outRow{col.count} = chanLoaded;
            outRow{col.file}  = fName;
            outRow{col.chan}  = unitNames{iu};
            
            % Get response property classification
            [pVals, pTypes, thisType, vmi, numSpks] = klGetType(spiketimes,Task);
            typeAlt = klGetTypeAlt(spiketimes,Task);
            resp = SPK_resp_profile(spiketimes,Task.SRT,1,0);
            [vPoiss, poissVisMat] = klPoissLatv2d(spiketimes,'-rwd',Task.Reward,'-stim',Task.StimOnset);
            [mPoiss, poissMovMat] = klPoissLatv2d(spiketimes-repmat(Task.SRT,1,size(spiketimes,2)),'-rwd',Task.Reward-Task.SRT,'-stim',Task.StimOnset-Task.SRT,'-stop',Task.SRT-Task.SRT,'direction','rev');
%             vPar = klTuneFit(vSDF(:,vTimes > 0 & vTimes < 200),Task.TargetLoc);
%             mPar = klTuneFit(mSDF(:,mTimes > -100 & mTimes < 0),Task.TargetLoc); 
            [vSig, vRise] = klSigLat(spiketimes);
            [mSig, mRise] = klSigLat(spiketimes-repmat(Task.SRT,1,size(spiketimes,2)),'-r',1,'base',[blMean,blStd]);
            if visMean < blMean, visSupp = 1; else visSupp = 0; end
            if movMean < blMean, movSupp = 1; else movSupp = 0; end
            if vSig > visMax, vSig = nan; vRise = nan; end
            if mSig < movMin, mSig = nan; mRise = nan; end
            
            % Cut visual latencies if it's a mov or a none, cut mov
            % latencies if it's a vis or a none
            if ismember(typeAlt,{'mov','none'}),
                [vPoiss, vSig, vRise] = deal(nan);
            elseif ismember(typeAlt,{'vis','none'}),
                [mPoiss, mSig, mRise] = deal(nan);
            end
            
%             vPar = klTuneEllipse(vSDF(:,vTimes > 0 & vTimes < 200),Task.TargetLoc);
%             mPar = klTuneEllipse(mSDF(:,mTimes > -100 & mTimes < 0),Task.TargetLoc); 
            
            % Get mean rates for target locations
            uLocs = unique(Task.TargetLoc); uLocs(isnan(uLocs)) = [];
            visLocRates = nan(1,length(uLocs));
            movLocRates = nan(1,length(uLocs));
            for il = 1:length(uLocs)
                visLocRates(il) = nanmean(nanmean(vSDF(Task.Correct == 1 & Task.TargetLoc == uLocs(il),vTimes >= vWind(1) & vTimes <= vWind(2)),2),1);
                movLocRates(il) = nanmean(nanmean(mSDF(Task.Correct == 1 & Task.TargetLoc == uLocs(il),mTimes <= -mWind(1) & mTimes >= -mWind(2)),2),1);
            end
            
            % Check if significantly visual
            try
                pVis = anovan(nanmean(vSDF(Task.Correct == 1,vTimes >= vWind(1) & vTimes <= vWind(2)),2),Task.TargetLoc(Task.Correct == 1),'display','off');
            catch pVis = nan;
            end
            try
                pMov = anovan(nanmean(mSDF(Task.Correct == 1,mTimes <= -mWind(1) & mTimes >= -mWind(2)),2),Task.TargetLoc(Task.Correct == 1),'display','off');
            catch pMov = nan;
            end
            
            if pVis < alph
                vPar = klTuneGauss(Task.TargetLoc(Task.Correct == 1),nanmean(vSDF(Task.Correct == 1,vTimes >= vWind(1) & vTimes <= vWind(2)),2),'-v',visualize);
            else
                vPar.mu = nan; vPar.sig = nan; vPar.amp = nan; vPar.bl = nan;
            end
            if pMov < alph,
                mPar = klTuneGauss(Task.TargetLoc(Task.Correct == 1),nanmean(mSDF(Task.Correct == 1,mTimes <= -mWind(1) & mTimes >= -mWind(2)),2),'-v',visualize);
            else
                mPar.mu = nan; mPar.sig = nan; mPar.amp = nan; mPar.bl = nan;
            end
            
            if ismember(typeAlt,{'mov','none'})
                vPar.mu = nan; vPar.sig = nan; vPar.amp = nan; vPar.bl = nan;
            end
            if ismember(typeAlt,{'vis','none'}),
                mPar.mu = nan; mPar.sig = nan; mPar.amp = nan; mPar.bl = nan;
            end
            
            if visualize,
                figure()
                subAx(1) = subplot(2,1,1);
                [sdf,sdfTimes]  = klSpkRatev2(spiketimes);
                pltMeanStd(sdfTimes(1):sdfTimes(2),nanmean(sdf,1),nanstd(sdf,[],1)./sqrt(size(sdf,1)),'k');
                vline(vPoiss);
                subAx(2) = subplot(2,1,2);
                for ir = 1:size(poissVisMat,1),
                    if sum(~isnan(poissVisMat(ir,:))) == 2,
                        plot(poissVisMat(ir,:),[ir,ir],'-k');
                        hold on;
                    end
                end
                linkaxes(subAx,'x');
%                 keyboard
            end
            
            outRow{col.type}  = thisType;
            outRow{col.type+1} = typeAlt{1};
            outRow{col.type+2} = typeAlt{2};
            outRow{col.type+3} = resp.resptype;
            outRow(col.mgPvals:(col.mgPvals+(length(pVals)-1))) = mat2cell(pVals,1,ones(1,length(pVals)));
            outRow{col.visIndex} = vmi(1);
            outRow{col.movIndex} = vmi(2);
            outRow{col.visPoiss}   = vPoiss;
            outRow{col.visPoiss+1} = (sum(~isnan(poissVisMat(Task.Correct,1)))/size(poissVisMat(Task.Correct,1),1))*100;
            outRow{col.movPoiss}   = mPoiss;
            outRow{col.movPoiss+1} = (sum(~isnan(poissMovMat(Task.Correct,1)))/size(poissMovMat(Task.Correct,1),1))*100;
            outRow{col.visLat}  = vSig;
            outRow{col.visLat+1} = vRise;
            outRow{col.movLat} = mSig;
            outRow{col.movLat+1} = mRise;
            outRow{col.visSupp}  = visSupp;
            outRow{col.movSupp}  = movSupp;
            outRow{col.vTune}    = vPar.mu;
            outRow{col.vTune+1}  = vPar.sig;
            outRow{col.vTune+2}  = vPar.amp;
            outRow{col.vTune+3}  = vPar.bl;
            outRow{col.mTune}    = mPar.mu;
            outRow{col.mTune+1}  = mPar.sig;
            outRow{col.mTune+2}  = mPar.amp;
            outRow{col.mTune+3}  = mPar.bl;
            
%             if mPar.sig > 1 || mPar.sig < 0, outRow{col.mTune+1} = nan; else outRow{col.mTune+1}  = mPar.sig; end
%             if vPar.sig > 1 || vPar.sig < 0, outRow{col.vTune+1} = nan; else outRow{col.vTune+1}  = vPar.sig; end
            
            
            % Check if the area string needs changed
            if sum(isletter(area)) == 0
                area = [' ',area];
            end
            area(strfind(area,'?')) = [];
            outRow{col.area} = area;
            chanNum     = str2double(DSPname(4:end-1));
            outRow{col.depthChan} = chanNum;
            
            %% Get spiking Statistics
            % {width,TfR,amplitude,meanRate,stdRate,fano,cv,cv2,lv,meanISI,stdISI,%ISI<3}
%            
%             wvMean = sortWaves.(monk{im}).(sprintf('s%s',sessNames{is}(~ismember(sessNames{is},'-')))).(unitNames{iu}).apWave;
%             wvTimes = sortWaves.(monk{im}).(sprintf('s%s',sessNames{is}(~ismember(sessNames{is},'-')))).(unitNames{iu}).apTimes;
            
            if isfield(meanSorts,monk{im}) && (myRow+4) <= length(meanSorts.(monk{im}).wave),
                wvMean = meanSorts.(monk{im}).wave{myRow+4};
                wvTimes = meanSorts.(monk{im}).time{myRow+4};
            else
                [wvMean, wvStd] = getMeanWaves(wave.waves,'-pm',0,'mean');
                wvTimes = ((1:32)-9).*25; 
            end
            
            % Get waveform width
%             [wvWidth, ~, ~, rpTime] = wv_width(wvMean);
            if size(wvMean,2) == 32, wvMeanUpSamp = 1; else wvMeanUpSamp = 0; end
            [wvWidth, rpTime] = klWvWidthv1(wvMean,wvTimes);%,'-t',wvTimes,'-u',wvMeanUpSamp);
            if isnan(wvWidth),
                wvSuccess = 0;
            else
                wvSuccess = 1;
            end
%             rpTime = nan;
            
            % Check for positive spikes
            isPos       = abs(max(wvMean)) > abs(min(wvMean));
            isPosThresh = wave.thresh > 0;
            maxAmp = max(abs(wvMean))*((-1)^isPos);
            relAmp = range(wvMean)./nanstd(wave.waves(:,1),1);
            
            spkShapeCell = {wvWidth,wvSuccess,rpTime,maxAmp,relAmp};
            outRow{col.SNR} = klGetSNRv1(wave.waves);
            
            %% Get rate statistics for all spikes
            [mnRate, stdRate] = klGetMeanRate(spiketimes);
            if blMean > rateLim,
                mnRate = nan; stdRate = nan; blMean = nan;
            end
            
%             fano              = SPK_get_Fano(spiketimes);
            fano              = klGetFano(spiketimes);
            cv                = klGetCV(spiketimes);
            cv2               = klGetCV(spiketimes,'-type','local');
            lv                = klGetLV(spiketimes);
            lvr               = klGetLV(spiketimes,'-type','revised');
            
            % Get ISI stats
%             blSpikes          = spiketimes; blSpikes(blSpikes > 0) = nan;
            isiMat            = diff(spiketimes,1,2);
            meanISI           = nanmean(isiMat(:));
            stdISI            = nanstd(isiMat(:));
            isiShort          = sum(isiMat(:) < 3)/sum(isfinite(isiMat(:)));
            
            spkCellT = {mnRate,stdRate,blMean,fano,cv,cv2,lv,lvr,meanISI,stdISI,isiShort};
            
            %% Now for baseline spikes
            blSpks = spiketimes;
            blSpks(blSpks < blWind(1) | blSpks > blWind(2)) = nan;
            
            [mnRate, stdRate] = klGetMeanRate(blSpks);
            if blMean > rateLim,
                mnRate = nan; stdRate = nan; blMean = nan;
            end
            
%             fano              = SPK_get_Fano(spiketimes);
            fano              = klGetFano(blSpks);
            cv                = klGetCV(blSpks);
            cv2               = klGetCV(blSpks,'-type','local');
            lv                = klGetLV(blSpks);
            lvr               = klGetLV(blSpks,'-type','revised');
            
            % Get ISI stats
%             blSpikes          = spiketimes; blSpikes(blSpikes > 0) = nan;
            isiMat            = diff(blSpks,1,2);
            meanISI           = nanmean(isiMat(:));
            stdISI            = nanstd(isiMat(:));
            isiShort          = sum(isiMat(:) < 3)/sum(isfinite(isiMat(:)));
            
            spkCellB = {mnRate,stdRate,blMean,fano,cv,cv2,lv,lvr,meanISI,stdISI,isiShort};
            
            % Put into a cell, then put that in outRow
            spkStatCell = {wvWidth,wvSuccess,rpTime,maxAmp,relAmp,mnRate,stdRate,blMean,fano,cv,cv2,lv,lvr,meanISI,stdISI,isiShort};
            outRow(col.spkShape:(col.spkShape+length(spkShapeCell)-1)) = spkShapeCell;
            outRow(col.spkStartT:(col.spkStartT+length(spkCellT)-1)) = spkCellT;
            outRow(col.spkStartB:(col.spkStartB+length(spkCellB)-1)) = spkCellB;
            
            outMat = cat(1,outMat,outRow);
            
            %[sdf, sdfTime] = klSpkRatev2(spiketimes);
            
        end
%         fprintf('Session %s analyzed in %s\n',sessNames{is},printTiming(sessTic));
%         keyboard
    end
    if ~isempty(outMat),
        xlswrite(xlFile,outMat,monk{im},sprintf('A%d',startRowOut+4'));
        fprintf('\nMonkey %s written in %s\n',monk{im},printTiming(monkTic));
    else
        fprintf('\nNo new data to write for Monkey %s\n',monk{im});
    end
%       keyboard
end

fprintf('All channels written in %s\n',printTiming(startTic));
