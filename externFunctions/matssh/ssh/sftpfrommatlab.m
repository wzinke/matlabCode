function status  =  sftpfrommatlab(userName,hostName,password,localfilename,remotefilename)
%sftpfrommatlab connects Matlab to a remote computer and uploads
%a file using the SFTP protocol
%
% STATUS  =  SFTPROMMATLAB(USERNAME,HOSTNAME,PASSWORD)
%
% Inputs:
%   USERNAME is the user name required for the remote machine
%   HOSTNAME is the name of the remote machine
%   PASSWORD is the password for the account USERNAME@HOSTNAME
%   LOCALFILENAME is the fully qualified path of the filename to be uploaded
%   REMOTEFILENAME is the fully qualified path where the file will be
%   stored at the remote computer
%
% Outputs:
%   STATUS is a Java com.trilead.ssh2.Connection object
%
% See also SSHFROMMATLAB,SSHFROMMATLABCLOSE, SSHFROMMATLABINSTALL, SSHFROMMATLABISSUE
%
% (c)2008 Athens Information Technology
%    Kostas Katrinis (kostas.katrinis@gmail.com)
%    Version 1.0
%

import com.trilead.ssh2.SFTPv3Client;
import com.trilead.ssh2.Connection;
import com.trilead.ssh2.Session;
import com.trilead.ssh2.SFTPv3FileHandle;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.BufferedInputStream;
import java.io.File;

%
%  Invocation checks
%
  if(nargin  ~=  5)
    error('Error: SftpROMMATLAB requires 5 input arguments...');
  end
  if(~ischar(userName)  || ~ischar(hostName)  ||  ~ischar(password) || ~ischar(localfilename) || ~ischar(remotefilename))
    error...
      (['Error: SftpROMMATLAB requires all input ',...
      'arguments to be strings...']);
  end


%Set up the connection with the remote server

try
        channel  =  Connection(hostName);
        channel.connect();
    catch
        error(['Error: SFTPFROMMATLAB could not connect to the'...
        ' remote machine %s ...'],...
        hostName);
end 

%
%  Check the authentication for login...
%
  
isAuthenticated = channel.authenticateWithPassword(userName,password);
if(~isAuthenticated)
    error...
        (['Error: SFTPFROMMATLAB could not authenticate the',...
        ' SSH connection...']);  
end

%Open session
channel2  =  channel.openSession();
sftpcl = SFTPv3Client(channel);

%create files
localf=sftpcl.createFile(localfilename);
remotefilename
remotef=sftpcl.openFileRW(remotefilename);

%transfer file byte by byte
  
br = BufferedInputStream(FileInputStream(localfilename));
count=0;
bufzi=0;

buf=zeros(1,1024);
bufsize=br.read(buf);
try
    while(bufsize~=-1)
        sftpcl.write(remotef,count,buf,0,bufsize);
        count=count+bufsize;
        bufsize=br.read(buf);
    end   
catch
    error(['Error: SFTPFROMMATLAB could not write to the'...
        ' remote machine %s ...'],...
        hostName);
end
count
sftpcl.close();
channel2.close();
channel.close();  

