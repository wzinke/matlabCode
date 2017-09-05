function sessSorts = tdtViewSorts(file,varargin)
clearvars -except file varargin

fresh = 1;

if ~exist('myDir','var'),
    myDir = 'Y:/Users/Kaleb/dataProcessed';
end
    
% Set constants
colors = 'rgbcmk';
sortStrs = {'PCA-KMeans','PCA-LSC','LPP-KMeans','LPP-LSC'};

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
    if exist(sprintf('%s/%s/Channel%d/autoSort_noAudit.mat',myDir,file,chans(ic)),'file'),
        if fresh || ~exist(sprintf('%s/%s/Channel%d/sortAudit1.mat',myDir,file,ic))
            % Pull this channel
    %         chanSorts = sessSortsTemp(chans(ic)).chanSorts;
            global subWaves chanSorts colors k mySet

            load(sprintf('%s/%s/Channel%d/autoSort_noAudit.mat',myDir,file,chans(ic)));
            
            % Extract info
            subWaves = chanSorts.alWaves(chanSorts.subInds,:);

            % Start sort loop
            for is = 1:4,
                figure(1);
                subplot(2,2,is);
                % Do unit loop
                for iu = 1:chanSorts.subSorts.allKs(is),
                    if is < 3,
                        scatter(chanSorts.pca(chanSorts.subSorts.idx{is}(:,chanSorts.subSorts.allKs(is))==iu,1),...
                            chanSorts.pca(chanSorts.subSorts.idx{is}(:,chanSorts.subSorts.allKs(is))==iu,2),[],colors(iu));
                        hold on;
                    else
                        scatter(chanSorts.lpp(chanSorts.subSorts.idx{is}(:,chanSorts.subSorts.allKs(is))==iu,1),...
                            chanSorts.lpp(chanSorts.subSorts.idx{is}(:,chanSorts.subSorts.allKs(is))==iu,2),[],colors(iu));
                        hold on;   
                    end
                end
                t=title(sortStrs{is},'fontsize',20);
                if is == chanSorts.subSorts.set,
                    set(t,'color',[.8 .2 .2]);
                end

                figure(10);
                subplot(2,2,is);

                % Do unit loop
                for iu = 1:chanSorts.subSorts.allKs(is),
                    plot(chanSorts.alTimes,subWaves(chanSorts.subSorts.idx{is}(:,chanSorts.subSorts.allKs(is))==iu,:),'color',colors(iu));
                    hold on;   
                end
                t=title(sortStrs{is});
                if is == chanSorts.subSorts.set,
                    set(t,'color',[.8 .2 .2]);
                end

            end

            good = 1;
            mySet = chanSorts.subSorts.set;
            k = chanSorts.subSorts.allKs(mySet);
            work = 0;
            sig = zeros(1,k);
            
            
            fprintf('To change set, type mySet = ...\n');
            fprintf('To change k, type k = ...\n');
            fprintf('If it needs more work, type good = 0\n');
            fprintf('To set a cluster as "unit", type sig(cluster) = 1\n');
            fprintf('\nTo continue: type work = 1, otherwise press F5\n');

            sig = sig(1:k);

            keyboard

            while work,
                figure(100);
                % Do unit loop
                for iu = 1:k,
                    if mySet < 3,
                        scatter(chanSorts.pca(chanSorts.subSorts.idx{mySet}(:,k)==iu,1),...
                            chanSorts.pca(chanSorts.subSorts.idx{mySet}(:,k)==iu,2),[],colors(iu));
                        hold on;
                    else
                        scatter(chanSorts.lpp(chanSorts.subSorts.idx{is}(:,k)==iu,1),...
                            chanSorts.lpp(chanSorts.subSorts.idx{is}(:,k)==iu,2),[],colors(iu));
                        hold on;   
                    end
                end
                t=title(sortStrs{mySet},'fontsize',20);
                set(t,'color',[.8 .2 .2]);

                figure(1000);
                % Do unit loop
                for iu = 1:k,
                    plot(chanSorts.alTimes,subWaves(chanSorts.subSorts.idx{mySet}(:,k)==iu,:),'color',colors(iu));
                    hold on;   
                end
                t=title(sortStrs{mySet});
                set(t,'color',[.8 .2 .2]);

                sig = zeros(1,k);

                fprintf('To change set, type mySet = ...\n');
                fprintf('To change k, type k = ...\n');
                fprintf('If no signal, type sig = 0\n');
                fprintf('If it needs more work, type good = 0\n');
                fprintf('\nTo accept: type work = 0, otherwise press F5\n');


                keyboard

                if k > length(sig),
                    sig = cat(2,sig,zeros(1,k-length(sig)));
                else
                    sig = sig(1:k);
                end
                close(100); close(1000);
            end

            sortAudit.set = mySet;
            sortAudit.k =k;
            sortAudit.good = good;
            sortAudit.isSig = sig;
            sortAudit.idx = chanSorts.subSorts.idx{mySet}(:,k);

            save(sprintf('%s/%s/Channel%d/sortAudit1.mat',myDir,file,ic),'sortAudit');

            sessSorts(chans(ic)) = sortAudit;
        else
            fprintf('Channel %d already sorted. Moving on\n',ic);
        end
            
    else
        fprintf('File %s, Channel %d not sorted. Moving on\n',file,ic);
    end
end

save(sprintf('%s/%s/sessSorts1.mat',myDir,file),'sessSorts');
   
end

function slowPlot()
    figure();
    global subWaves chanSorts colors k mySet
    for i = 1:k,
        plot(subWaves(chanSorts.subSorts.idx{mySet}(:,k)==i,:)','color',colors(i)); hold on; pause;
    end
end
            
            
            
            