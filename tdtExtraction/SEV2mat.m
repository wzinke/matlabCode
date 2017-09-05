function [data] = SEV2mat(SEV_DIR, varargin)
%SEV2MAT  TDT SEV file format extraction.
%   data = SEV2mat(SEV_DIR), where SEV_DIR is a string, retrieves
%   all sev data from specified directory in struct format. SEV files
%   are generated by an RS4 Data Streamer, or by setting the Unique
%   Channel Files option in Stream_Store_MC or Stream_Store_MC2 macro
%   to Yes.
%
%   data    contains all continuous data (sampling rate and raw data)
%
%   data = SEV2mat(SEV_DIR,'parameter',value,...)
%
%   'parameter', value pairs
%      'CHANNEL'    integer, returns the sev data from specified channel
%                       only (default = 0 for all channels)
%      'JUSTNAMES'  boolean, retrieve only the valid event names
%      'EVENTNAME'  string, specific event name to retrieve data from
%      'VERBOSE'    boolean, set to false to disable console output
%      'DEVICE'     string, connect to specific RS4 device.  DEVICE can be
%                       the IP address or NetBIOS name of RS4-device
%                       (e.g. RS4-41001).  Requires TANK and BLOCK
%                       parameters
%      'TANK'       string, tank on RS4 to retrieve data from. Requires
%                       DEVICE and BLOCK parameters
%      'BLOCK'      string, block on RS4 to retrieve data from. Requires
%                       DEVICE and TANK parameters

if ~mod(nargin, 2)
    error('not enough input arguments')
end

% defaults
CHANNEL   = 0;
EVENTNAME = '';
DEVICE    = '';
TANK      = '';
BLOCK     = '';
T1        = 0;
T2        = 0;
VERBOSE   = 1;
JUSTNAMES = 0;

% parse varargin
for i = 1:2:length(varargin)
    eval([upper(varargin{i}) '=varargin{i+1};']);
end

if any([~isempty(DEVICE) ~isempty(TANK) ~isempty(BLOCK)])
    if any([isempty(DEVICE) isempty(TANK) isempty(BLOCK)])
        error('DEVICE, TANK and BLOCK must all be specified');
    else
        SEV_DIR = sprintf('\\\\%s\\data\\%s\\%s\\', DEVICE, TANK, BLOCK);
    end
end

data = [];

ALLOWED_FORMATS = {'single','int32','int16','int8','double','int64'};

if strcmp(SEV_DIR(end), filesep) == 0
    SEV_DIR = [SEV_DIR filesep];
end

file_list = dir([SEV_DIR '*.sev']);
nfiles = length(file_list);
if nfiles < 1
    warning(['no sev files found in ' SEV_DIR])
    return
end

% find out what data we think is here
for i = 1:length(file_list)
    [pathstr, name, ext] = fileparts(file_list(i).name);
    
    % find channel number
    matches = regexp(name, '_[Cc]h[0-9]*', 'match');
    if ~isempty(matches)
        sss = matches{end};
        file_list(i).chan = str2double(sss(4:end));
    end
    
    % find starting hour
    matches = regexp(name, '-[0-9]*h', 'match');
    if ~isempty(matches)
        sss = matches{end};
        file_list(i).hour = str2double(sss(2:end-1));
    else
        file_list(i).hour = 0;
    end
    
    % check file size
    file_list(i).data_size = file_list(i).bytes - 40;
    
    path = [SEV_DIR file_list(i).name];
    fid = fopen(path, 'rb');
    if fid < 0
        warning([path ' not opened'])
        return
    end
    
    % create and fill streamHeader struct
    streamHeader = [];
    
    streamHeader.fileSizeBytes   = fread(fid,1,'uint64');
    streamHeader.fileType        = char(fread(fid,3,'char')');
    streamHeader.fileVersion     = fread(fid,1,'char');
    
    % event name of stream
    s = regexp(file_list(i).name, '_', 'split');
    streamHeader.eventName = s{end-1};
    
    if streamHeader.fileVersion < 3
        
        %if streamHeader.fileVersion == 2
        oldEventName  = char(fread(fid,4,'char')');
        %else
        %    streamHeader.eventName  = fliplr(char(fread(fid,4,'char')'));
        %end
        
        % current channel of stream
        streamHeader.channelNum        = fread(fid, 1, 'uint16');
        % total number of channels in the stream
        streamHeader.totalNumChannels  = fread(fid, 1, 'uint16');
        % number of bytes per sample
        streamHeader.sampleWidthBytes  = fread(fid, 1, 'uint16');
        reserved                 = fread(fid, 1, 'uint16');
        
        % data format of stream in lower four bits
        streamHeader.dForm      = ALLOWED_FORMATS{bitand(fread(fid, 1, 'uint8'),7)+1};
        
        % used to compute actual sampling rate
        streamHeader.decimate   = fread(fid, 1, 'uint8');
        streamHeader.rate       = fread(fid, 1, 'uint16');
    else
        error(['unknown version ' num2str(streamHeader.fileVersion)]);
    end
    
    % compute sampling rate
    if streamHeader.fileVersion > 0
        streamHeader.fs = 2^(streamHeader.rate)*25000000/2^12/streamHeader.decimate;
    else
        % make some assumptions if we don't have a real header
        streamHeader.dForm = 'single';
        streamHeader.fs = 24414.0625;
        warning('%s has empty header; assuming ch %d format %s and fs = %.2f\nupgrade to OpenEx v2.18 or above\n', ...
            file_list(i).name,  ...
            file_list(i).chan, streamHeader.dForm, streamHeader.fs);
    end
    
    % check variable name
    %varname = matlab.lang.makeValidName(streamHeader.eventName);
    varname = streamHeader.eventName;
    for ii = 1:numel(varname)
        if ii == 1
            if ~isnan(str2double(varname(ii)))
                varname(ii) = 'x';
            end
        end
        if ~isletter(varname(ii)) && isnan(str2double(varname(ii)))
            varname(ii) = '_';
        end
    end
    
    if ~isvarname(streamHeader.eventName)
        warning('%s is not a valid Matlab variable name, changing to %s', streamHeader.eventName, varname);
    end
    
    func = str2func(streamHeader.dForm);
    tempvar = func(zeros(1,1));
    w = whos('tempvar');
    file_list(i).npts = file_list(i).data_size / w.bytes;
    file_list(i).fs = streamHeader.fs;
    file_list(i).dForm = streamHeader.dForm;
    file_list(i).eventName = streamHeader.eventName;
    file_list(i).varName = varname;
    fclose(fid);
end

eventNames = unique({file_list.eventName});
if JUSTNAMES
    data = eventNames;
    return
end

for ev = 1:numel(eventNames)
    
    thisEvent = eventNames{ev};
    
    if ~strcmp(EVENTNAME, '') && ~strcmp(EVENTNAME, thisEvent)
        continue
    end
    
    file_list_temp = [];
    for j = 1:length(file_list)
        if strcmp(file_list(j).eventName, thisEvent)
            file_list_temp = [file_list_temp file_list(j)];
        end
    end
    
    max_chan = max([file_list_temp.chan]);
    max_hour = max([file_list_temp.hour]);
    hour_values = sort(unique([file_list_temp.hour]));
    
    % preallocate data array
    if CHANNEL > 0
        matching_ch = find([file_list_temp.chan] == CHANNEL);
    else
        matching_ch = find([file_list_temp.chan] == 1);
    end
    
    % determine offsets if there is chunking
    total_samples = 0;
    ct = 1;
    offsets = ones(1, numel(hour_values));
    for jjj = hour_values
        temp_num = intersect(find([file_list_temp.hour] == jjj), matching_ch);
        total_samples = total_samples + file_list_temp(temp_num).npts;
        ct = ct + 1;
        offsets(ct) = total_samples;
    end
    
    % now allocate it
    if CHANNEL > 0
        data.(file_list_temp(1).varName).data(1, total_samples) = func(0);
        loop = CHANNEL;
    else
        data.(file_list_temp(1).varName).data(max_chan, total_samples) = func(0);
        loop = 1:max_chan;
    end
    
    % loop through the channels
    
    for ii = 1:length(loop)
        i = loop(ii);
        matching_ch = find([file_list_temp.chan] == i);
        
        % loop through the chunks
        for j = hour_values
            
            file_num = intersect(find([file_list_temp.hour] == j), matching_ch);
            
            % open file
            path = [SEV_DIR file_list_temp(file_num).name];
            fid = fopen(path, 'rb');
            if fid < 0
                warning([path ' not opened'])
                return
            end
            
            % skip first 40 bytes
            fread(fid, 10, 'single');
            
            % read rest of file into data array as correct format
            varname = file_list_temp(file_num).varName;
            
            data.(varname).name = streamHeader.eventName;
            data.(varname).fs = streamHeader.fs;
            
            this_offset = offsets(find(hour_values == j));
            data.(varname).data(ii,this_offset:this_offset+file_list_temp(file_num).npts-1) = fread(fid, inf, ['*' streamHeader.dForm])';
            
            % close file
            fclose(fid);
            
            if VERBOSE
                file_list(file_num)
            end
        end
    end
end
