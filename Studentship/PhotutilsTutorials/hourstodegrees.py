#converts HH:MM:SS format to decimal format for degrees
#also antiquated

def hoursToDeg(instring:str):
    list1 = instring.split(" ")
    RA = list1[0]
    Dec = list1[1]
    list2 = RA.split(':')
    list3 = Dec.split(':')
    hours = int(list2[0])
    minutes = int(list2[1])
    seconds = float(list2[2])
    if seconds > 0:
        totalSeconds = seconds + 60*minutes + 60*60*hours
    else:
        totalSeconds = seconds - 60*minutes - 60*60*hours
    RADecimal = totalSeconds * 360/86400
    degrees = int(list3[0])
    minutes2 = int(list3[1])
    seconds2 = float(list3[2])
    if degrees > 0:
        totalSeconds2 = 60*minutes2+seconds2
    else:
        totalSeconds2 = -1*(60*minutes2+seconds2)
    DECDecimal = degrees + totalSeconds2/3600
    return (RADecimal, DECDecimal)

print(hoursToDeg('06:29:08.430 -71:29:38.98'))