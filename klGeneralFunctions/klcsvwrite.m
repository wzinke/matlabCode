function success = klcsvwrite(fileName,in,varargin)

% Set defaults
append = 0;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-h'}
            header = varargin{varStrInd(iv)+1};
        case {'-a'}
            append = varargin{varStrInd(iv)+1};
    end
end

success = 0;
try
    % Open the file and generate fid
    fid = fopen(fileName,'w+');

    % Write the header if it exists
    if exist('header','var')
        for ic = 1:(length(header)-1)
            switch isnumeric(header{ic})
                case 0
                    fprintf(fid,'%s,',header{ic});
                case 1
                    fprintf(fid,'%.5f,',header{ic});
            end
        end
        switch isnumeric(header{end})
            case 0
                fprintf(fid,'%s\n',header{end});
            case 1
                fprintf(fid,'%.5f\n',header{end});
        end
    end
    
    % Write the data
    for ir = 1:size(in,1)
        for ic = 1:(size(in,2)-1)
            switch isnumeric(in{ir,ic})
                case 0
                    fprintf(fid,'%s,',in{ir,ic});
                case 1
                    fprintf(fid,'%.5f,',in{ir,ic});
            end
        end
        switch isnumeric(in{ir,end})
            case 0
                fprintf(fid,'%s\n',in{ir,end});
            case 1
                fprintf(fid,'%.5f\n',in{ir,end});
        end
    end
    
    success = 1;
end


% Close fid
fclose(fid);