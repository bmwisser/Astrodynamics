def TimeToUTC(year, month, day, hour, min, sec, UTC_offset):
    import math
    """Takes a given date and time and adjusts it to UTC time, making corrections to the date if necessary."""
    #Check that the UTC Offset does not exceed 24 hours, this is non-sensical
    if abs(UTC_offset) > 24:
        print("Error: UTC offset specification is greater than 24 hours. Please re-adjust to a sensical value. "
              "Exiting...")
        return None
    def month_limit(month, year = 0):
        import math
        if month == 2 and year == 0:
            print("ERROR! If the month is February, the year must also be input. Exiting...")
            return None
        #Determine the day limit of the month based on the input year and month
        if month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12:
            day_limit = 31
        elif month == 4 or month == 6 or month == 9 or month == 11:
            day_limit = 30
        elif month == 2:
            remainder = year % 4
            if remainder == 0:
                day_limit = 29
            else:
                day_limit = 28
        return day_limit

    #Check if UTC offset input is an integer or some fractional hour
    frac_check = UTC_offset - math.floor(UTC_offset)
    if frac_check != 0: #Hour fraction! Must add some minutes
        min_add = frac_check * 60 #convert from hours to min
        min = min + math.floor(min_add)
        frac_check = min_add - math.floor(min_add)
        if frac_check != 0: #Minute fraction! Must add some seconds
            sec_add = frac_check * 60 #convert from min to seconds
            sec = sec + sec_add

    hour = hour + math.floor(UTC_offset)
    #If the hour exceeds the bounds for a day, adjust the day and hour numbers
    if hour > 24:
        day = day + 1
        hour = hour % 24
    if hour < 0:
        day = day - 1
        hour = 24 + hour
    #Check Day Bounds for the month
    day_limit = month_limit(month, year)
    if day < 1:
        month = month - 1
        if month < 1:
            month = 12
            year = year - 1
            day = month_limit(month, year)
        else:
            day = month_limit(month, year)
    elif day > day_limit:
        month = month + 1
        if month > 12:
            month = 1
            year = year + 1
        day = 1

    return year, month, day, hour, min, sec
def JulianDate(year,month,day,hour,min,sec,UTC_offset=0):
    import math
    """Calculates the Julian Date of a given epoch."""
    #Check that the inputs fit within the valid bounds
    if year >= 2099 or year <=1901:
        print("Error: Year exceeds acceptable range of 1901-2099. Exiting...")
        return None
    if month < 1 or month > 12:
        print("Error: Month falls outside of acceptable range 1-12. Exiting...")
        return None
    if day < 1 or day > 31:
        print("Error: Day falls outside of acceptable range 1-31. Exiting...")
        return None

    #Adjust the input time to UTC time. Capture any potential day changes that can occur.
    if UTC_offset != 0:
        year, month, day, hour, min, sec = TimeToUTC(year, month, day, hour, min, sec, UTC_offset)
    J0 = 367*year - math.floor(7/4*(year+math.floor((month+9)/12)))+math.floor(275*month/9)+day+1721013.5
    DayFrac = (hour + min/60 + sec/3600) / 24

    JulianDate = J0 + DayFrac
    return JulianDate, J0


# for dd in [18,22,23]:
#     print(dd)
#     [JD,J0] = JulianDate(2018,10,dd,21,00,00,6)
#     print(JD)