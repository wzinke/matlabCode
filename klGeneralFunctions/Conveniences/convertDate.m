%% Takes date and switches format between m/d/y and yymmdd

function [dateOut thisYear thisMonth thisDay] = convertDate(dateIn)

switch class(dateIn)
    
    case 'double'
        dateAsString = num2str(dateIn);
        if length(dateAsString) > 6, fprintf('\nDate %s is too long...\n\tReturning to command line\n',dateAsString); return; end;
        thisYear = dateAsString(1:2);
        thisMonth = dateAsString(3:4); if strcmpi(thisMonth(1),'0'), thisMonth = thisMonth(2); end;
        thisDay = dateAsString(5:6); if strcmpi(thisDay(1),'0'), thisDay = thisDay(2); end;
        
        dateOut = sprintf('%s/%s/%s',thisMonth,thisDay,thisYear);
       
    case 'char'
        slashInd = strfind(dateIn,'/'); if length(slashInd) > 2, fprintf('Date is in improper format...\n\tReturning to command line\n'); return; end;
        thisMonth = dateIn(1:slashInd(1)-1); if length(thisMonth) < 2, thisMonth = sprintf('0%s',thisMonth); end;
        thisDay = dateIn(slashInd(1)+1:slashInd(2)-1); if length(thisDay) < 2, thisDay = sprintf('0%s',thisDay); end;
        thisYear = dateIn(slashInd(2)+1:end); if length(thisYear) > 2, thisYear = thisYear(end-1:end); end;
        
        dateOutStr = sprintf('%s%s%s',thisYear,thisMonth,thisDay);
        dateOut = str2num(dateOutStr);
end