function newDate = subtractDate(startDate,dayIncrement)

if dayIncrement >= 0
    
    [dateAsNum dateYear dateMonth dateDay] = convertDate(startDate);
    daysPerMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
    
    newDay = str2num(dateDay)-dayIncrement;
    newMonth = str2num(dateMonth);
    newYear = str2num(dateYear);
    
    daysLastMonth = daysPerMonth(newMonth-1);
    daysIncrMonth = daysLastMonth;
    
    while newDay <= 0
        newDay = daysIncrMonth+newDay;
        newMonth = newMonth-1;
        if newMonth <= 0, newYear = newYear - 1; newMonth = length(daysPerMonth); end;
        if newMonth == 1
            daysIncrMonth = daysPerMonth(length(daysPerMonth));
        else
            daysIncrMonth = daysPerMonth(newMonth-1);
        end
    end
    
    newYearStr = num2str(newYear); 
    newMonthStr = num2str(newMonth); if length(newMonthStr) < 2, newMonthStr = sprintf('0%s',newMonthStr); end;
    newDayStr = num2str(newDay); if length(newDayStr) < 2, newDayStr = sprintf('0%s',newDayStr); end;
    newDateStr = sprintf('%s%s%s',newYearStr,newMonthStr,newDayStr);

    newDate = str2num(newDateStr);

else
    newDate = addDate(startDate,-dayIncrement);
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    