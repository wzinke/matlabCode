function [lat, nb]= poissLat_Longv2(spks,Task)

visualize = 1;

%% Get the time ranges for each trial (using -500 as the minimum time)
tRange = Task.Reward - (-500);
spks(spks > repmat(Task.Reward,1,size(spks,2))) = nan;
numSpks = sum(isfinite(spks),2);
r  = numSpks./tRange;
isi = diff(spks,1,2);
meanISI = nanmean(isi,2);

%% Trial by trial
for it = 1:size(spks,1),
    numBursts = 1;
    
    thisTrSpks = spks(it,:);
    thisTrISI  = isi(it,:);
    thisISIthresh = meanISI(it)/2;
    
    % thisISIthresh = meanISI/2 as per Legendy & Salcman 1985
    isiUnderThresh = thisTrISI < thisISIthresh;
    startLook = find(thisTrSpks > Task.StimOnset(it),1);
    keepLooking = 1;
    while keepLooking
        foundBurst = 0;
        for ii = startLook:(length(isiUnderThresh)-2)
            % Find the first set of 3 spikes with mean ISI / thresh
            if length(unique(isiUnderThresh(ii:ii+2))) == 1 && unique(isiUnderThresh(ii:ii+2)) == 1,
                foundBurst = 1;
                bStartTent = ii;
                % Get surprise index for this set of spikes
                startP = poisspdf(2,r(it)*(spks(it,ii+2)-spks(it,ii)));
                startS = -log(startP);
                startSpk = ii;
                break
            end
        end
        if foundBurst == 0,
            keepLooking = 0;
            break
        end
        % Set starting vars
        oldS = startS;
        newS = startS;
        newSpks = 2;
        keepGoing = 1;
        while keepGoing == 1,
            % Does including an extra spike increase surprise?
            while newS >= oldS && thisTrISI(startSpk+newSpks) <= 2*meanISI(it),
                oldS = newS;
                newSpks = newSpks+1;
                newP = poisspdf(newSpks,r(it)*(spks(it,startSpk+newSpks)-spks(it,startSpk)));
                newS = -log(newP);
            end
            keepGoing = 0;
            tempEnd   = newSpks;
            % Above loop stops when newS stops increasing surprise.
            % Does this hold true for (up to) 10 more spikes?
            compS = oldS;
            extraSpks = 0;
            while (extraSpks < 10 && (thisTrISI(startSpk+newSpks) <= 2*meanISI(it))) && ((newSpks + startSpk) <= length(thisTrISI))
                if (startSpk + newSpks) == 21,
                    %keyboard
                end
                newP = poisspdf(newSpks,r(it)*(spks(it,startSpk+newSpks)-spks(it,startSpk)));
                newS = -log(newP);
                if newS > compS, 
                    keepGoing = 1; 
                    break; 
                end;
                newSpks = newSpks + 1;
                extraSpks = extraSpks + 1;
            end
        end

        % Set new start vars
        endSpk = startSpk+newSpks;
        latestS = newS;
        subSpks = 0;
        checkSpks = endSpk-startSpk;

        % Does removing a spike from the beginning of the burst increase
        % surprise?
        while newS >= latestS && subSpks < checkSpks,
            subSpks = subSpks + 1;
            newP = poisspdf(checkSpks-subSpks,r(it)*(spks(it,endSpk) - spks(it,startSpk+subSpks)));
            newS = -log(newP);
        end

        finalStart = startSpk + (subSpks-1);
        finalEnd   = endSpk;
        
        finalP = poisspdf(sum(thisTrSpks >= thisTrSpks(finalStart) & thisTrSpks >= thisTrSpks(finalEnd)),r(it)*(thisTrSpks(finalEnd) - thisTrSpks(finalStart)));
        if finalP < .05,
            burst{it,numBursts} = [spks(it,finalStart), spks(it,finalEnd)];
        end
        
        numBursts = numBursts + 1;
        
        startLook = finalEnd;
        %keyboard

        if visualize && finalP < .05
            figure(); hold on;
            scatter(spks(it,:),repmat(it,size(spks,2),1),'k*');
            plot([burst{it,numBursts-1}(1), burst{it,numBursts-1}(2)],[it+2 it+2],'r');
            %plot([bStart, bEnd],[it+4, it+4],'b');
            legend('Spks','Activity','Burst');
            set(gca,'YLim',[it - 10, it + 10]);
            t = title(sprintf('%d Burst %d - Final p = %.3f, si = %.3f',it,numBursts-1,finalP,-log(finalP)));
            keyboard
        end
         %keyboard
    end
    %keyboard
end
%keyboard
[lat{1:size(burst,2)}] = deal(cell(1,size(burst,2)));
nb = nan(1,size(burst,2));
for ib = 1:size(burst,2)
    mnLat{ib} = nanmean(cell2mat(burst(:,ib)));
    medLat{ib} = nanmedian(cell2mat(burst(:,ib)));
    nb(ib)  = sum(~cellfun(@isempty,burst(:,ib)));
end
keyboard