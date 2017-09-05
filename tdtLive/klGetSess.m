function Task = klGetSess(sessName,varargin)

% Set defaults
print = 0;
nPrint = 500;
fresh = 0;
% rawDir = '/mnt/teba/data/Kaleb/antiSessions/';
rawDir = '/mnt/teba/Users/Kaleb/proAntiRaw';
procDir = '/mnt/teba/users/Kaleb/proAntiProcessed/';

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-p','print'}
            print = varargin{varStrInd(iv)+1};
        case {'-f','fresh'}
            fresh = varargin{varStrInd(iv)+1};
    end
end


if ~ismember(sessName(end),['/','\']), sessName(end+1) = filesep;  end
if ismember(rawDir(end),['/','\']), fileName = [rawDir,sessName]; else fileName = [rawDir,filesep,sessName]; end
if ismember(procDir(end),['/','\']), saveDir = [procDir,sessName]; else saveDir = [procDir,filesep,sessName]; end

if ~exist(saveDir,'file')
    mkdir(saveDir);
end

if ispc
    getFun = @TDT2mat;
elseif isunix
    getFun = @TDTbin2mat;
end

% Pull out event codes and times
if exist(sprintf('%sevsRaw.mat',saveDir),'file') && ~fresh
   load(sprintf('%sevsRaw.mat',saveDir));
else
    if print
        printStr = 'Getting Events...';
        fprintf(printStr);
    end

    tdtEvRaw = getFun(fileName,'TYPE',{'epocs','scalars'},'VERBOSE',0);
    if isfield(tdtEvRaw.epocs,'STRB')
        tdtEvs = tdtEvRaw.epocs.STRB.data;
        if any(tdtEvs > 2^15)
            tdtEvs = tdtEvs-2^15;
        end
        tdtEvTms = tdtEvRaw.epocs.STRB.onset.*1000;
        tdtEvTms(tdtEvs <= 0) = [];
        tdtEvs(tdtEvs <= 0) = [];
    else
    %     tdtEvRaw = getFun(fileName,'TYPE',{'scalars'},'VERBOSE',0);
        if isempty(tdtEvRaw.scalars)
            Task = [];
            return
        end
        tdtEvs = tdtEvRaw.scalars.EVNT.data;
        if any(tdtEvs >= (2^15))
            tdtEvs = tdt2EvShft(tdtEvs);
        end
        if any(mod(tdtEvs,1)) ~= 0
            tdtEvs = tdtEvRaw.scalars.EVNT.data - (2^15);
        end
        tdtEvTms = tdtEvRaw.scalars.EVNT.ts.*1000;
        tdtEvTms(tdtEvs < 0) = [];
        tdtEvs(tdtEvs < 0) = [];
    end
    save(sprintf('%sevsRaw.mat',fileName),'tdtEvs','tdtEvTms');
end

% Get eye positions
if exist(sprintf('%sEyes.mat',saveDir),'file') && ~fresh
    load(sprintf('%sEyes.mat',saveDir));
else
%     try
        if print
            for ib = 1:length(printStr)
                fprintf('\b');
            end
            printStr = 'Getting Eye X...';
            fprintf(printStr);
        end
        eyeXRaw = getFun(fileName,'TYPE',{'streams'},'STORE','EyeX','VERBOSE',0);
        if print
            for ib = 1:length(printStr)
                fprintf('\b');
            end
            printStr = 'Getting Eye Y...';
            fprintf(printStr);
        end
        eyeYRaw = getFun(fileName,'TYPE',{'streams'},'STORE','EyeY','VERBOSE',0);
        eyeX = eyeYRaw.streams.EyeY.data.*(3);
        eyeY = eyeXRaw.streams.EyeX.data.*(-3);
        eyeT = (0:(length(eyeX)-1)).*(1000/eyeXRaw.streams.EyeX.fs);
        eyeR = sqrt(eyeX.^2 + eyeY.^2);
        eyeThRaw = klRad2Deg(atan(abs(eyeY)./abs(eyeX)));
        eyeTh = nan(size(eyeThRaw));
        eyeTh(eyeX > 0 & eyeY > 0) = eyeThRaw(eyeX > 0 & eyeY > 0);
        eyeTh(eyeX < 0 & eyeY > 0) = 180-eyeThRaw(eyeX < 0 & eyeY > 0);
        eyeTh(eyeX < 0 & eyeY < 0) = 180+eyeThRaw(eyeX < 0 & eyeY < 0);
        eyeTh(eyeX > 0 & eyeY < 0) = 360-eyeThRaw(eyeX > 0 & eyeY < 0);
        Eyes.X = eyeX;
        Eyes.Y = eyeY;
        Eyes.R = eyeR;
        Eyes.Theta = eyeTh;
        Eyes.Times = eyeT;
        if print
            for ib = 1:length(printStr)
                fprintf('\b');
            end
            printStr = 'Saving Eyes...';
            fprintf(printStr);
            save(sprintf('%sEyes.mat',saveDir),'Eyes');
        end
        clear eyeXRaw eyeYRaw eyeX eyeY eyeT eyeR eyeThRaw eyeTh
        Eyes.Good = 1;
%     catch
%         Eyes.X = nan;
%         Eyes.Y = nan;
%         Eyes.R = nan;
%         Eyes.Theta = nan;
%         Eyes.Times = nan;
%         Eyes.Good = 0;
%     end
end
Task = klGetTask(tdtEvs,tdtEvTms,Eyes);
if print
    fprintf('\n');
    printStr = 'Saving Behavior...';
    fprintf(printStr);
end
save(sprintf('%sBehav.mat',saveDir),'Task');
delete(sprintf('%sevsRaw.mat',saveDir));