function status  =  scpfrommatlab(userName,copyFromhostName,copyTohostName,password,copyFromfilename,copyTofilename)
%sftpfrommatlab connects Matlab to a remote computer and uploads
%a file using SCP. Assume existing SCP implementation installed in
%the local system that is called from within Matlab. (e.g. Putty in this
%case)
%
% STATUS  =  SCPFROMMATLAB(USERNAME,copyFromhostName,copyTohostName,PASSWORD,copyFromfilename,copyTofilename)
%
% Inputs:
%   USERNAME is the user name required for the remote machine
%   copyFromhostName is the name of machine where the file that will be
%   transfered from resides. Use 'localhost' to imply that the local
%   system if the source host.
%   copyTohostName is the name of machine where the file will be transfered
%   to. Use 'localhost' to imply that the local
%   system if the destination host.
%   PASSWORD is the password for the account USERNAME@HOSTNAME
%   copyFromfilename is the fully qualified path of the filename to be uploaded
%   copyTofilename is the fully qualified path where the file will be
%   stored at the remote computer
%
% Outputs:
%   STATUS: 0 if transfer completed successfully.
%           1 if transfer failed.
%
% (c)2008 Athens Information Technology
%    Kostas Katrinis (kkat@ait.edu.gr)
%    Version 1.0
%

% specify name of the scp executable in the local system
% NOTE: the full path where the scp executable resides
%       should have been added to the Matlab path.
load('sshfrommatlab2Configure.mat');

%
%  Invocation checks
%
  if(nargin  ~=  6)
    error('Error: SCPFROMMATLAB requires 6 input arguments...');
  end
  if(~ischar(userName)  || ~ischar(copyFromhostName) || ~ischar(copyTohostName) ||  ~ischar(password) || ~ischar(copyFromfilename) || ~ischar(copyTofilename))
    error...
      (['Error: SCPFROMMATLAB requires all input ',...
      'arguments to be strings...']);
  end

%invoke local implementation of SCP
if strcmp(copyFromhostName,'localhost')==0
    if strcmp(copyTohostName,'localhost')==0
        command = sprintf('%s -q -l %s -pw %s %s:%s %s:%s',scpLocalPath,userName,password,copyFromhostName,copyFromfilename,copyTohostName,copyTofilename);
    else
        command = sprintf('%s -q -l %s -pw %s %s:%s %s',scpLocalPath,userName,password,copyFromhostName,copyFromfilename,copyTofilename);
    end
else
    command = sprintf('%s -q -l %s -pw %s %s %s:%s',scpLocalPath,userName,password,copyFromfilename,copyTohostName,copyTofilename);
end
[status,result] = system(command);
