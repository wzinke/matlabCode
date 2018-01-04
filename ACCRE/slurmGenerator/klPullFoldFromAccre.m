function ranThrough = klPullFoldFromAccre(tarName,userOpts,varargin),

% Set defaults
remoteDir = '.';
localDir = '.';
localTars = './tempTars';
print = 1;

% Decode varargin
varStrInd=find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-r'},
            remoteDir = varargin{varStrInd(iv)+1};
        case {'-l'},
            localDir = varargin{varStrInd(iv)+1};
        case {'-t'},
            tarDir = varargin{varStrInd(iv)+1};
    end
end

if ~exist('tarDir','var'),
    tarDir = remoteDir;
end

if exist('userInfo.mat') && ischar(userOpts),
    load('userInfo.mat');
    tmpStruct.name = userOpts;
    tmpStruct.pass = passes{ismember(users,userOpts)};
    userOpts = tmpStruct; clear tmpStruct;
end

ranThrough = 0;
if ~exist(localDir,'file'),
    mkdir(localTars);
end

% try
    % Open connection
    if print, fprintf('Opening Connection...\n'); end
    sshStruct = ssh2_config('login.accre.vanderbilt.edu',userOpts.name,userOpts.pass);

    if print, fprintf('Making .tar file...\n'); end
    sshStruct = ssh2_command(sshStruct, sprintf('tar -C %s -cf %s/%s.tar %s .',remoteDir,tarDir,tarName,tarName));
    
    % Send the folder
    if print, fprintf('Copying folder...\n'); end
    sshStruct = scp_get(sshStruct,sprintf('%s.tar',tarName),localTars,tarDir);
    
    % Close connection temporarily...
    sshStruct = ssh2_close(sshStruct);
    
    % Unzip it locally
    if print, fprintf('Unzipping .tar file...\n'); end
    untar(sprintf('%s/%s.tar',localTars,tarName),sprintf('%s',localDir));

    % Delete the remote .tar to save space
    if print, fprintf('Removing unnecessary files...\n'); end
    sshStruct = ssh2_config('login.accre.vanderbilt.edu',userOpts.name,userOpts.pass);
    sshStruct = ssh2_command(sshStruct, sprintf('rm %s/%s.tar',tarDir,tarName));
    delete(sprintf('%s/%s.tar',localTars,tarName));
    
    % Close connection
    sshStruct = ssh2_close(sshStruct);
    
    if print, fprintf('Copying successful!\n'); end
    ranThrough = 1;
% end
       