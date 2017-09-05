%% klBookKeepingNew


% Initialize some variables
monk = {'Gauss','Helmholtz'};
xlFile = './klDataBookKeeping_mg.xlsx';
dataFold = 'Y:/Users/Wolf/ephys_db';

type        = 'DSP';                             % DSP for spikes, LFP for LFP, and EEG for EEG
task        = 'MG';
physType    = 'rates';

rOff = 4;
cOff = 3;

startCell = [num2abc(cOff+1),num2str(rOff+1)];
load xlCols

for im = 1:length(monk)
    thisMonk =monk{im};
    
    % Initialize Variables
    outMat = {};
    
    %% Start Loops
    allSess = dir(sprintf('%s/%s/2*',dataFold,thisMonk));
    sessNames = {allSess.name};
    for is = 1:length(sessNames),
        sessUnits = dir(sprintf('%s/%s/%s/%s/*DSP*',dataFold,thisMonk,sessNames{is},type));
        unitNames = {sessUnits.name};
        
        chanLoaded = 0;
        for iu = 1:length(unitNames),
            clear f fName
            outRow = cell(1,40);
            
            % Check for valid data files
            f = dir(sprintf('%s/%s/%s/%s/%s/%s_%s_%s*',dataFold,monk{im},sessNames{is},type,unitNames{iu},sessNames{is},unitNames{iu},task));
            if length(f) == 0,
                fprintf('\tNo "%s" file found for unit %s... Skipping to next unit\n',task,unitNames{iu});
                continue;
            elseif length(f) > 1,
                fprintf('\tMultiple "%s" files found for unit %s... Please check directory structure\n',task,unitNames{iu});
                %keyboard
            end
            fName = f(1).name;
                
            % Load unit data
            loadTic = tic;
            
            load(sprintf('%s/%s/%s/%s/%s/%s',dataFold,monk{im},sessNames{is},type,unitNames{iu},fName));
            load(sprintf('%s/%s/%s/%s/%s/waves/%s_waves.mat',dataFold,monk{im},sessNames{is},type,unitNames{iu},fName(1:(end-4))));
            
            %% Get Identification Info
            outRow{col.sess}  =  sessNames{is};
            chanLoaded = chanLoaded + 1;
            outRow{col.count} = chanLoaded;
            outRow{col.file}  = fName;
            outRow{col.chan}  = unitNames{iu};
            
            % Get response property classification
            [pVals, pTypes, thisType, vmi, numSpks] = klGetType(spiketimes,Task);
            resp = SPK_resp_profile(spiketimes,Task.SaccEnd);
            
            outRow{col.type}  = thisType;
            outRow(col.mgPvals:(col.mgPvals+(length(pVals)-1))) = mat2cell(pVals,1,ones(1,length(pVals)));
            outRow{col.visIndex} = vmi(1);
            outRow{col.movIndex} = vmi(2);
            
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
            [wvMean, wvStd] = getMeanWaves(wave.waves,'-pm',0,'median');
            
            % Get waveform width
            [wvWidth, ~, ~, rpTime] = wv_width(wvMean,'-t',1:length(wvMean));
            
            % Check for positive spikes
            maxAmp = max([abs(wvMean)]);
            isPos       = abs(max(wvMean)) > abs(min(wvMean));
            isPosThresh = wave.thresh > 0;
            
            % Get rate statistics
            [mnRate, stdRate] = klGetMeanRate(spiketimes);
            fano              = SPK_get_Fano(spiketimes);
            cv                = klGetCV(spiketimes);
            cv2               = klGetCV(spiketimes,'-type','local');
            lv                = klGetLV(spiketimes);
            
            % Get ISI stats
            isiMat            = diff(spiketimes,1,2);
            meanISI           = nanmean(isiMat(:));
            stdISI            = nanstd(isiMat(:));
            isiShort          = sum(isiMat(:) < 3)/sum(isfinite(isiMat(:)));
            
            % Put into a cell, then put that in outRow
            spkStatCell = {wvWidth,rpTime,maxAmp,mnRate,stdRate,fano,cv,cv2,lv,meanISI,stdISI,isiShort};
            outRow(col.spkStart:(col.spkStart+length(spkStatCell)-1)) = spkStatCell;
            outMat = cat(1,outMat,outRow);
            
            %[sdf, sdfTime] = klSpkRatev2(spiketimes);
            
        end
        %keyboard
    end
    xlswrite(xlFile,outMat,monk{im},'A5');
end
        