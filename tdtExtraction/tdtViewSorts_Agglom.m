function sessSorts = tdtViewSorts_Agglom(file,varargin)
clearvars -except file varargin

fresh = 1;

if ~exist('myDir','var'),
    myDir = 'Y:/Users/Kaleb/dataProcessed';
end
    
% Set constants
colors = 'rgbcmk';
myPols = {'neg','pos'};

% Get number of channels for loop
chanDir = dir(sprintf('%s/%s/Channel*',myDir,file));
chanNames = {chanDir.name};
nChans = length(chanNames);
chans = nan(1,nChans);
for i = 1:nChans,
    chans(i) = str2num(chanNames{i}(8:end));
end
chans = sort(chans);

% Load task if necessary
if nChans > 0,
    load(sprintf('%s/%s/Behav.mat',myDir,file));
end

% Initialize output
sessSorts(max(chans)) = struct('set',[],'k',[],'good',[],'isSig',[],'idx',[]);

% Load temporary sorts
% x=tic;
% load(sprintf('%s/sessSortsTemp.mat',file));
% fprintf('sessSortsTemp loaded in %s\n',printTiming(x));
% global sessSortsTemp
% 
% Start channel loop
for ic = 1:nChans,
    fprintf('Showing channel %d...\n',chans(ic));
    chanFold = sprintf('%s/%s/Channel%d',myDir,file,chans(ic));
    close all;
    clear sortAudit chanSorts subWaves
    if exist(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',myDir,file,chans(ic)),'file'),
        % Pull this channel
%         chanSorts = sessSortsTemp(chans(ic)).chanSorts;
        global subWaves chanSorts colors k idx

        load(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',myDir,file,chans(ic)));

        % Extract info
        for ip = 1:2,
            close all
            subWaves = chanSorts.(myPols{ip}).allWaves(chanSorts.(myPols{ip}).subInds,:);
            
            % Plot waveforms
            figure(1);
            for ik = 1:chanSorts.(myPols{ip}).k,
                plot(subWaves(chanSorts.(myPols{ip}).idx==ik,:)','color',colors(ik)); hold on;
            end
            
            % Plot scatterplot
            figure(2);
            for ik = 1:chanSorts.(myPols{ip}).k,
                scatter3(chanSorts.(myPols{ip}).subScores(chanSorts.(myPols{ip}).idx==ik,1),chanSorts.(myPols{ip}).subScores(chanSorts.(myPols{ip}).idx==ik,2),chanSorts.(myPols{ip}).allTimes(chanSorts.(myPols{ip}).subInds(chanSorts.(myPols{ip}).idx==ik)),[],colors(ik)); hold on;
            end
            
            work = 0;
            k = chanSorts.(myPols{ip}).k;
            idx = chanSorts.(myPols{ip}).idx;
            sig = zeros(1,k);
            
            fprintf('To change k, type k = ...\n');
            fprintf('If it needs more work, type good = 0\n');
            fprintf('To set a cluster as "unit", type sig(cluster) = 1\n');
            fprintf('\nTo continue: type work = 1, otherwise press F5\n');

            keyboard
            
            while work,
                idx=map2idx(chanSorts.(myPols{ip}).aggMap,k);
                figure(100);
                for ik = 1:k,
                    plot(subWaves(idx==ik,:)','color',colors(ik)); hold on;
                end
                
                figure(101);
                for ik = 1:k,
                    scatter3(chanSorts.(myPols{ip}).subScores(idx==ik,1),chanSorts.(myPols{ip}).subScores(idx==ik,2),chanSorts.(myPols{ip}).allTimes(chanSorts.(myPols{ip}).subInds(idx==ik)),[],colors(ik)); hold on;
                end
                
                fprintf('To change k, type k = ...\n');
                fprintf('If it needs more work, type good = 0\n');
                fprintf('To set a cluster as "unit", type sig(cluster) = 1\n');
                fprintf('\nTo continue: type work = 1, otherwise press F5\n');

                keyboard
                
                if k > length(sig),
                    sig = cat(2,sig,zeros(1,k-length(sig)));
                else
                    sig = sig(1:k);
                end
                
                close(100); close(101);
                
            end       
            
            chanSorts.(myPols{ip}).k = k;
            chanSorts.(myPols{ip}).idx = idx;
            chanSorts.(myPols{ip}).sig = sig;
                
        end
        
        save(sprintf('%s/%s/Channel%d/autoSortAgglom_noAudit.mat',myDir,file,chans(ic)),'chanSorts','-v7.3');
    else
        fprintf('File %s, Channel %d not sorted. Moving on\n',file,ic);
    end
end
   
end

function slowPlot()
    figure();
    global subWaves colors k idx
    for i = 1:k,
        plot(subWaves(idx==i,:)','color',colors(i)); hold on; pause;
    end
end
            
function newIDX = map2idx(inMap,k)
    % Now make an easy to access summary/companion matrices
    outCount = nan(size(inMap));
    outInds = nan(size(inMap));
    for i = 1:length(inMap),
        uClusts = unique(inMap(:,i));
        [outCount(1:length(uClusts),i), sortInds] = sort(hist(inMap(:,i),uClusts),'descend');
        outInds(1:length(uClusts),i) = uClusts(sortInds);
    end
        
    myCol = find(outCount(k,:) == max(outCount(k,:)),1,'last');
    uK = outInds(1:k,myCol);
    newIDX = nan(size(inMap,1),1);
    for ik = 1:length(uK),
        newIDX(inMap(:,myCol)==uK(ik)) = ik;
    end    
end
            