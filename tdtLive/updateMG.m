%%% Update MG
function updateMG(lookChans,tdt)

close all;

if nargin < 2
    tdt = SynapseLive();
end

% Hold on to old data
oldData = tdt.data;
% Get new data
tdt.update;
% Concatenate data structure
tdt.data = klCatStruct(oldData,tdt.data);

% Get event code information
[tdtEvs, tdtEvTms] = checkEvs(tdt.data);
Task = klGetTask(tdtEvs,tdtEvTms);

close all;
for ic = 1:length(lookChans)
    fprintf('Getting MG Responses for Channel %d (%d of %d): ',lookChans(ic),ic,length(lookChans));
    chanSpikes = tdt.data.snips.eNe1.ts(tdt.data.snips.eNe1.chan==lookChans(ic)).*1000;
    
    klMGPhysio(Task,chanSpikes,'-c',lookChans(ic));
    fprintf('Done...\n');
%     figure(1); figure(2);
    
%     keyboard
    
%     close all;
end

end


function [tdtEvs, tdtEvTms] = checkEvs(dataRaw)
    if isfield(dataRaw.epocs,'STRB')
        tdtEvs = dataRaw.epocs.STRB.data;
        if any(tdtEvs > 2^15)
            tdtEvs = tdtEvs-2^15;
        end
        tdtEvTms = dataRaw.epocs.STRB.onset.*1000;
        tdtEvTms(tdtEvs <= 0) = [];
        tdtEvs(tdtEvs <= 0) = [];
    else
    %     tdtEvRaw = getFun(fileName,'TYPE',{'scalars'},'VERBOSE',0);
        if isempty(dataRaw.scalars)
            tdtEvs = nan;
            tdtEvTms = nan;
            return
        end
        tdtEvs = dataRaw.scalars.EVNT.data;
        if any(tdtEvs >= (2^15))
            tdtEvs = tdt2EvShft(tdtEvs);
        end
        if any(mod(tdtEvs,1)) ~= 0
            tdtEvs = dataRaw.scalars.EVNT.data - (2^15);
        end
        tdtEvTms = dataRaw.scalars.EVNT.ts.*1000;
        tdtEvTms(tdtEvs < 0) = [];
        tdtEvs(tdtEvs < 0) = [];
    end
end