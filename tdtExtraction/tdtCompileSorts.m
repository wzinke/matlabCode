function sessSortsTemp = tdtCompileSorts(file)

% Get number of channels for loop
chanDir = dir(sprintf('%s/Channel*',file));
chanNames = {chanDir.name};
nChans = length(chanNames);
chans = nan(1,nChans);
for i = 1:nChans,
    chans(i) = str2num(chanNames{i}(8:end));
end
chans = sort(chans);

% Initialize structure
sessSortsTemp(max(chans)) = struct('chanSorts',[]);

% Start channel loop
for ic = 1:nChans,
    fprintf('Fetching channel %d...\n',chans(ic));
    chanFold = sprintf('%s/Channel%d',file,chans(ic));
    close all;
    clear sortAudit chanSorts subWaves
    if exist(sprintf('%s/autoSort_noAudit.mat',chanFold),'file'),
        load(sprintf('%s/autoSort_noAudit.mat',chanFold));

        sessSortsTemp(chans(ic)).chanSorts = chanSorts;
    end
end

save(sprintf('%s/sessSortsTemp.mat',file),'sessSortsTemp');
