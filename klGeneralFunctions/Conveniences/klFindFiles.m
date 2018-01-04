function [fileFolds,fileNames] = klFindFiles(root,ext)

if ~strcmpi(ext(1),'.')
    ext = ['.',ext];
end

% First, check this base folder
checkTop = dir([root,'/*',ext]);
fileFolds(1:length(checkTop),1) = {checkTop.folder};
fileNames(1:length(checkTop),1) = {checkTop.name};

% Now, check for subfolders
subCheck = dir(root);
subFolds = {subCheck.name};
subFolds = subFolds([subCheck.isdir]);
subFolds = subFolds(cellfun(@(x) ~strcmpi(x(1),'.'),subFolds));

for is = 1:length(subFolds)
    [subFold,subNames] = klFindFiles([root,'/',subFolds{is}],ext);
    fileFolds = cat(1,fileFolds,subFold);
    fileNames = cat(1,fileNames,subNames);
end