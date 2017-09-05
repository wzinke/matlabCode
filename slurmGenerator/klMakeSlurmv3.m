function [slurmFileDir, slurmFileName] = klMakeSlurmv3(mFile,varargin)

%% Set defaults
slurmFileName = 'matlab.slurm';
nCores = 1;
matVers = 'r2016a';
memAmnt = 500;
memPerCore = ceil(memAmnt./nCores);
memUnits = 'M';
outFile = 'matlab_job_slurm.out';
nNodes = 1;
timeStr = '00:10:00';
nTasks = 1;
slurmFileDir = '.';
doGit = 0;
gitLink = '/users/placeholder/alsoaplaceholder';
mailUser = 'kaleb.a.lowe@vanderbilt.edu';
mailType = 'ALL';
asFun = 0;
% mTargDir = '';

%% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)},
        case {'-c','ncores','nCores'},
            nCores = varargin{varStrInd(iv)+1};
        case {'-s','sFile'},
            slurmFileName = varargin{varStrInd(iv)+1};
        case {'sDir'},
            slurmFileDir = varargin{varStrInd(iv)+1};
        case {'-g','gitLink','gitlink'},
            gitLink = varargin{varStrInd(iv)+1};
            doGit = 1;
        case {'-m'},
            memAmnt = varargin{varStrInd(iv)+1};
		case {'-mpc'},
			memPerCore = varargin{varStrInd(iv)+1};
        case {'-f'},
            asFun = 1;
            funIn = varargin{varStrInd(iv)+1};
        case {'-t'},
            timeStr = varargin{varStrInd(iv)+1};
		case {'-mu'},
			memUnits = varargin{varStrInd(iv)+1};
		case {'-out'},
			outFile = varargin{varStrInd(iv)+1};
		case {'-mt'},
			mailType = varargin{varStrInd(iv)+1};
    end
end

%% Save SLURM parameters as a structure in order to load it into a script
slurmParams.cores = nCores;
slurmParams.timeStr = timeStr;
slurmParams.version = matVers;
slurmParams.memory.amount = memAmnt;
slurmParams.memory.units = memUnits;
slurmParams.nTasks = nTasks;
slurmParams.nNodes = nNodes;
slurmDots = strfind(slurmFileName,'.');
save(sprintf('%s/slurmParams_%s.mat',slurmFileDir,slurmFileName(1:(slurmDots(end)-1))),'slurmParams','-v7.3');

%% Make sure mFile has .m suffix
dotInd = strfind(mFile,'.');
if isempty(dotInd), dotInd = length(mFile)+1; end

% Cut last period to end of file name, then attach '.m'
mFile = [mFile(1:(dotInd(end)-1)),'.m'];

%% First, open the file
sFile = fopen(sprintf('%s/%s',slurmFileDir,slurmFileName),'w+');

%% Now print the slurm commmands: (needs further editing/flexibility
fprintf(sFile,'#!/bin/bash\n');
fprintf(sFile,sprintf('#SBATCH --mail-user=%s\n',mailUser));
fprintf(sFile,sprintf('#SBATCH --mail-type=%s\n',mailType));
% fprintf(sFile,sprintf('#SBATCH --nodes=%d\n',nNodes));
% fprintf(sFile,sprintf('#SBATCH --tasks-per-node=%d\n',nCores));
% fprintf(sFile,sprintf('#SBATCH --ntasks=%d\n',nTasks));
fprintf(sFile,sprintf('#SBATCH --time=%s\n',timeStr));
% if nCores > 1,
% 	fprintf(sFile,sprintf('#SBATCH --mem-per-core=%d%s\n',memPerCore,memUnits));
% else
	fprintf(sFile,sprintf('#SBATCH --mem=%d%s\n',memAmnt,memUnits));
% end
fprintf(sFile,sprintf('#SBATCH --output=%s\n',outFile));
fprintf(sFile,'\n');
if doGit,
    fprintf(sFile,sprintf('setpkgs -a git\n'));
    fprintf(sFile,sprintf('git clone %s\n',gitLink));
end
fprintf(sFile,sprintf('setpkgs -a matlab_%s\n',matVers));
fprintf(sFile,'\n');

if asFun,
    fprintf(sFile,sprintf('matlab -nodisplay -nosplash -r "addpath(genpath(''.'')),%s(%s),quit()"',mFile(1:(end-2)),funIn));
else
    fprintf(sFile,sprintf('matlab -nodisplay -nosplash < %s\n',mFile));
end

fclose(sFile);

