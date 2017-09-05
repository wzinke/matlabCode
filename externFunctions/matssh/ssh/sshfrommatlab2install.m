function sshfrommatlabinstall
%SSHFROMMATLAB2INSTALL provides installation instructions for SSHFROMMATLAB
%
% Installation Instructions
% -------------------------
%
% The Matlab code included in this distribution is a frontend for calling
% SSH/SCP/SFTP functions, whose actual protocol implementation is
% implemented by external tools. 
%
% Specifically, you will need following implementations installed in your
% system (freeware at the time of writing):
%
% 1. SSH and SFTP implementation in Java called Trilead SSH for Java. To install
% Trilead SSH for Java:
%
%   a) Copy the respective jar file into a directory under your Matlab root, 
%      e.g. under C:\MATLAB7\java\jar\com\
%   b) Append the fully qualified path (including file name) of the jar file
%   into your Matlab's "classpath.txt" file, which is usually found under
%   "/your Matlab root folder/toolbox/local".
%   (Example entry under Windows XP:
%   C:\MATLAB7\java\jar\com\trilead-ssh2-build213.jar)
%
% 2. If you also need SCP functionality, you will need an SCP implementation 
% installed in the native system, e.g. PuTTy in
% MS Window OS. Install the implemented software and make sure that you add
% environment variables to your OS, such that scp and sftp are called with
% the shortcut commands from your command prompt/shell. Note that you do
% not need to follow steps 2 and 3, if you can't do without SCP, i.e. if
% you only need SSH and/or SFTP.
%
% 3. Open the .mat file 'sshfrommatlab2Configure.mat' and set the
% scpLocalPath value to the fully qualified command (or just command if it
% is already in your path variables) name of the native systems SCP
% executable. E.g. with PuTTY on MS Windows XP and path variables set, this
% translates to setting:
% scpLocalPath=pscp
%
% (c)2008 Athens Information Technology
%    Kostas Katrinis (kkat@ait.edu.gr)
%    Version 1.0
%
help sshfrommatlab2install;