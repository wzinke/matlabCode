rawDir = 'Y:/Users/Kaleb/dataRaw';
procDir = 'Y:/Users/Kaleb/dataProcessed';
inFile = 'Init_SetUp-160711-151215';
load(sprintf('%s/%s/manChanThreshes.mat',procDir,inFile));
    
for ic = 1:32,
    % Read in data
    chanSEVs = SEV2mat_kl(sprintf('%s/%s',rawDir,inFile),'CHANNEL',ic,'VERBOSE',0);
    
    if ~any(chanSEVs.Wav1.data > 1),
        chanSEVs.Wav1.data = chanSEVs.Wav1.data.*1000000;
        chanThresh(ic) = chanThresh(ic).*1000000;
    end
    
    figure();
    plot(chanSEVs.Wav1.data(1:100000));
    hline(chanThresh(ic));
    
    keyboard
    close all
end