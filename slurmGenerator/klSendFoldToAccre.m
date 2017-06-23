function ranThrough = klSendFoldToAccre(foldPath,tarName,userOpts,varargin),

% Set defaults
remoteDir = '.';
localDir = '.';
print = 1;

% Decode varargin
varStrInd=find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd),
    switch varargin{varStrInd(iv)},
        case {'-rdir'},
            remoteDir = varargin{varStrInd(iv)+1};
        case {'-ldir'},
            localDir = varargin{varStrInd(iv)+1};
    end
end

ranThrough = 0;
if ~exist(localDir,'file'),
    mkdir(localDir);
end

% try
    % Make the .tar file
    if print, fprintf('Making .tar file...\n'); end
    tar(sprintf('%s/%s.tar',localDir,tarName),foldPath);

    % Open connection
    if print, fprintf('Opening Connection...\n'); end
    sshStruct = ssh2_config('login.accre.vanderbilt.edu',userOpts.name,userOpts.pass);

    % Send the folder
    if print, fprintf('Copying folder...\n'); end
    sshStruct = ssh2_command(sshStruct, sprintf('mkdir %s',remoteDir));
    sshStruct = scp_put(sshStruct,sprintf('%s.tar',tarName),remoteDir,localDir,sprintf('%s.tar',tarName));

    % Unzip it remotely
    if print, fprintf('Unzipping...\n'); end
%     sshStruct = ssh2_command(sshStruct, sprintf('mkdir ./%s',tarName));
    sshStruct = ssh2_command(sshStruct, sprintf('tar -C %s -xf %s/%s.tar',remoteDir,remoteDir,tarName));

    % Delete the original .tar in order to save some space
    if print, fprintf('Removing unnecessary files...\n'); end
    sshStruct = ssh2_command(sshStruct, sprintf('rm %s/%s.tar',remoteDir,tarName));

    % Close connection
    sshStruct = ssh2_close(sshStruct);

    % Delete local tar folder
    delete(sprintf('%s/%s.tar',localDir,tarName));
    
    if print, fprintf('Copying successful!\n'); end
    ranThrough = 1;
% end
       