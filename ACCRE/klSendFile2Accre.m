function success = klSendFile2Accre(fileDir,varargin)

% Set defaults
remoteRoot = '/home/loweka';
direction = 1; % 1 = send, 2 = pull

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-u','user','-user'},
            accreUser = varargin{varStrInd(iv)+1};
        case {'-p','pass','password'},
            accrePass = varargin{varStrInd(iv)+1};
        case {'-d'},
            direction = varargin{varStrInd(iv)+1};
        case {'-r'},
            remoteRoot = varargin{varStrInd(iv)+1};
    end
end

success = 0;
try
    % so "fileDir" is the local, let's make a string for the remote
    [path,name,ext] = fileparts(fileDir);
    remoteStr = sprintf('%s@login.accre.vanderbilt.edu:%s/%s',accreUser,remoteRoot,name);

    % Check to see if this is a directory
    if isempty(ext)
        isDir = 1;
    else
        isDir = 0;
        remoteStr = [remoteStr,ext];
    end

    % Let's make the sshpass statement header
    statement = sprintf('!sshpass -p ''%s'' scp ',accrePass);

    % Let's add the recursion if we need to:
    if isDir
        statement = [statement, '-r '];
    end

    % Now finish the statement based on directionality
    if direction == 1,
        statement = sprintf('%s -o StrictHostKeyChecking=no %s %s ',statement,fileDir,remoteStr);
    else
        statement = sprintf('%s -o StrictHostKeyChecking=no %s %s',statement,remoteStr,fileDir);
    end

    eval(statement);
    success = 1;
end