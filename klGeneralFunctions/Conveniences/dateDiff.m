function diff = dateDiff(day1,day2)

if day1 < day2
    largeDay = day2;
    smallDay = day1;
elseif day1 >= day2
    largeDay = day1;
    smallDay = day2;
end


checkDay = smallDay;
daysAdded = -1;
while checkDay < largeDay
    daysAdded = daysAdded + 1;
    checkDay = addDate(smallDay,daysAdded);
end

diff = daysAdded;