function klMakeSlurmv1(mFile,varargin)

%% Set defaults
slurmFileName = 'matlab.slurm';
nCores = 1;
matVers = 'r2014a';
memAmnt = 500;
memUnits = 'M';
outFile = 'matlab_job_slurm.out';
nNodes = 1;
timeStr = '00:10:00';
nTasks = 1;
slurmFileDir = '.';

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
    end
end

%% Make sure mFile has .m suffix
dotInd = strfind(mFile,'.');
if isempty(dotInd), dotInd = length(mFile)+1; end

% Cut last period to end of file name, then attach '.m'
mFile = [mFile(1:(dotInd(end)-1)),'.m'];

%% First, open the file
sFile = fopen(sprintf('%s/%s',slurmFileDir,slurmFileName),'w+');

%% Now print the slurm commmands: (needs further editing/flexibility
fprintf(sFile,'#!/bin/bash\n');
fprintf(sFile,sprintf('#SBATCH --nodes=%d\n',nNodes));
fprintf(sFile,sprintf('#SBATCH --ntasks=%d\n',nTasks));
fprintf(sFile,sprintf('#SBATCH --time=%s\n',timeStr));
fprintf(sFile,sprintf('#SBATCH --mem=%d%s\n',memAmnt,memUnits));
fprintf(sFile,sprintf('#SBATCH --output=%s\n',outFile));
fprintf(sFile,'\n');
fprintf(sFile,sprintf('setpkgs -a matlab_%s\n',matVers));
fprintf(sFile,'\n');
fprintf(sFile,sprintf('matlab -nodisplay -nosplash < %s\n',mFile));

fclose(sFile);

