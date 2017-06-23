thisFile = 'Init_SetUp-160829-153146';
for ic = 1:64,
	load(strcat(rootDir,'/',thisFile,'/Channel',num2str(ic),'/autoSortAgglom_noAudit.mat'));
	chanSorts.neg=tdtSortWavesv1a(chanSorts.neg.allTimes,chanSorts.neg.allWaves,'-d',3);
	chanSorts.pos=tdtSortWavesv1a(chanSorts.pos.allTimes,chanSorts.pos.allWaves,'-d',3);
	save(strcat(rootDir,'/',thisFile,'/Channel',num2str(ic),'/autoSortAgglom_noAudit.mat'));
end
