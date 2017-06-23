function chanSorts = klChanSortv1(thisFile,ic,varargin)

% Set defaults
rootDir = '.';
maxWvs = 40000;
nDims = 2;

% Decode varargin
varStrInd=find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-r'},
            rootDir = varargin{varStrInd(iv)+1};
        case {'-m','-w'},
            maxWvs = varargin{varStrInd(iv)+1};
        case {'-d','nDims'},
            nDims = varargin{varStrInd(iv)+1};
    end
end
addpath(genpath(rootDir));

fprintf('\n\n\nmaxWvs = %d before loading...\n\n',maxWvs);


load(strcat(rootDir,'/',thisFile,'/Channel',num2str(ic),'/autoSortAgglom_noAudit.mat'),'-mat','chanSorts');
z=fieldnames(chanSorts.neg);
fprintf('chanSorts.neg fields:\n');
for iz = 1:length(z),
	fprintf('\t%s\n',z{iz});
end
chanSorts.neg=tdtSortWavesv2a(chanSorts.neg.allTimes,chanSorts.neg.allWaves,'-d',nDims,'-w',maxWvs,'-s',['_Ch',num2str(ic)],'tempFold','/scratch/loweka');
chanSorts.pos=tdtSortWavesv2a(chanSorts.pos.allTimes,chanSorts.pos.allWaves,'-d',nDims,'-w',maxWvs,'-s',['_Ch',num2str(ic)],'tempFold','/scratch/loweka');
save(strcat(rootDir,'/',thisFile,'/Channel',num2str(ic),'/autoSortAgglom_noAudit.mat'),'chanSorts','-v7.3');
