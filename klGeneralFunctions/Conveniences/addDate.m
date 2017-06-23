function newDate = addDate(startDate,dayIncrement)

%keyboard;
if dayIncrement >= 0
    [dateAsNum dateYear dateMonth dateDay]= convertDate(startDate);
    daysPerMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
    newDay = str2num(dateDay) + dayIncrement;
    newMonth = str2num(dateMonth);
    newYear = str2num(dateYear);

    daysThisMonth = daysPerMonth(newMonth);
    daysIncrMonth = daysThisMonth;
    while newDay > daysIncrMonth
        newDay = newDay-daysIncrMonth;
        newMonth = newMonth+1;
        if newMonth > 12, newMonth = newMonth - 12; newYear = newYear + 1; end
        daysIncrMonth = daysPerMonth(newMonth);
    end

    %if newDay > daysThisMonth, newDay = mod(newDay,daysThisMonth); newMonth = newMonth + 1; end;
    %newYear = str2num(dateYear); if newMonth > 12, newMonth = 1; newYear = newYear + 1; end;

    newYearStr = num2str(newYear); 
    newMonthStr = num2str(newMonth); if length(newMonthStr) < 2, newMonthStr = sprintf('0%s',newMonthStr); end;
    newDayStr = num2str(newDay); if length(newDayStr) < 2, newDayStr = sprintf('0%s',newDayStr); end;
    newDateStr = sprintf('%s%s%s',newYearStr,newMonthStr,newDayStr);

    newDate = str2num(newDateStr);

else
    newDate = subtractDate(startDate,-dayIncrement);
end
%keyboard
%blockCheckCutoff = sprintf('%s%s%s',dateYear,num2str(
