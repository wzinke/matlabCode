function mFileOut = tdtMakeSessSlurmv1(thisFile,nDims)

if nargin < 2,
    nDims = 2;
end

mFileName = thisFile;
mFileName(strfind(mFileName,'-')) = '_';
mFile = fopen(sprintf('./klSortSlurm_%s.m',mFileName),'w+');
mFileOut = sprintf('./klSortSlurm_%s.m',mFileName);

fprintf(mFile,sprintf('thisFile = ''%s'';\n',thisFile));
fprintf(mFile,sprintf('for ic = 1:64,\n'));
fprintf(mFile,'\tload(strcat(rootDir,''/'',thisFile,''/Channel'',num2str(ic),''/autoSortAgglom_noAudit.mat''));\n');
fprintf(mFile,sprintf('\tchanSorts.neg=tdtSortWavesv1a(chanSorts.neg.allTimes,chanSorts.neg.allWaves,''-d'',%d);\n',nDims));
fprintf(mFile,sprintf('\tchanSorts.pos=tdtSortWavesv1a(chanSorts.pos.allTimes,chanSorts.pos.allWaves,''-d'',%d);\n',nDims));
fprintf(mFile,'\tsave(strcat(rootDir,''/'',thisFile,''/Channel'',num2str(ic),''/autoSortAgglom_noAudit.mat''));\n');
fprintf(mFile,'end\n');

fclose(mFile);