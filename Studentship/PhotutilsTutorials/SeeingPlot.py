from random import gauss
from tokenize import Single
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.io import fits
from projectMultiple import imageProject
from stackMultiple import imageStack
from FinalAperturePhotometry import* 
from psfForStackedImages import*
from astropy.stats import sigma_clip, sigma_clipped_stats, gaussian_sigma_to_fwhm
import os
from psfForStackedImages import SingleImagePSF
import datetime
import matplotlib.dates as mdates
import math
pixToArc = 0.8305


def plotSeeing(objects:str, dates:str, method:str, mode):

    """
    objects is a list of objects separated by commas. JL163,WD0830,JL163,WD1153.....\n
    can have repeated objects\n
    dates - 20221214,20230112,20230115.... are the corresponding dates to those objects\n
    method - mean or median\n
    mode - 1 = short timescale, unstacked images over a single night. Uses reprojected images. Will be used with one single date\n
         - 2 = longer timescale, stacked images (8 for whole night) and compare different nights\n
    
    """
    dateList = dates.split(',')
    objectList = objects.split(',')
    
    if mode == 1:
        MedianSigmaB = []
        MedianSigmaV = []
        BTimeList = []
        VTimeList = []
        BYearList = []
        BMonthList = []
        BDayList = []
        VYearList = []
        VMonthList = []
        VDayList = []
        i = 0
        for item in objectList:
            fileB = f'Studentship/PhotutilsTutorials/90secExposures{item}_B_{dateList[i]}.txt'
            if os.path.exists(f'Studentship/PhotutilsTutorials/90secExposures{item}_V_{dateList[i]}.txt'):
                fileV = f'Studentship/PhotutilsTutorials/90secExposures{item}_V_{dateList[i]}.txt'
            else:
                fileV = fileB
    
            with open(fileB, 'r') as f1:
                lines = f1.readlines()
            for line in lines:
                line = line.strip('\n')
                header = fits.getheader(line, 0)
                startTime:str = header['DATE-OBS']
                timeString = (startTime.strip('\'').split('T'))
                time = timeString[1]
                timeSplit = time.split(':')
                hours = int(timeSplit[0])
                minutes = int(timeSplit[1])
                seconds = int(math.floor(float(timeSplit[2])))
                realHours = (int(hours + 8)) % 24
                dateSplit = (timeString[0]).split('-')
                day = int(dateSplit[2])
                
                if (hours >= 16) and (minutes>0):
                    day = day + 1
                    print(day)
                BYearList.append(int(dateSplit[0]))
                BMonthList.append(int(dateSplit[1]))
                BDayList.append(day)         
            
                BTimeList.append(datetime.datetime(int(dateSplit[0]), int(dateSplit[1]), day, realHours, minutes, seconds))
                print(line)
                sigmaList = SingleImagePSF(line, item)
                MedianSigmaB.append(np.median(sigmaList))
            
            with open(fileV, 'r') as f2:
                lines2 = f2.readlines()
            for line2 in lines2:
                line2 = line2.strip('\n')
                header2 = fits.getheader(line2, 0)
                startTime2:str = header2['DATE-OBS']
                timeString2 = (startTime2.strip('\'').split('T'))
                time2 = timeString2[1]
                timeSplit2 = time2.split(':')
                hours2 = int(timeSplit2[0])
                minutes2 = int(timeSplit2[1])
                seconds2 = int(math.floor(float(timeSplit2[2])))
                realHours2 = (int(hours2 + 8)) % 24
                dateSplit = (timeString2[0]).split('-')
                day = int(dateSplit[2])
                
                if (hours2 >= 16) and (minutes2 >= 0):
                    day = day + 1
                    print(day)
                VYearList.append(int(dateSplit[0]))
                VMonthList.append(int(dateSplit[1]))
                VDayList.append(day)
                
                VTimeList.append(datetime.datetime(int(dateSplit[0]), int(dateSplit[1]), day, realHours2, minutes2, seconds2))
                sigmaList2 = SingleImagePSF(line2, item)
                MedianSigmaV.append(np.median(sigmaList2))
            i += 1

        Bdatetimes = BTimeList
        Vdatetimes = VTimeList
        #i = 0
        #ii = 0
        #for t in BTimeList:
            #Bdatetimes.append(datetime.datetime.combine(datetime.date(BYearList[i], BMonthList[i], BDayList[i]), t))
            #i += 1
        #for j in VTimeList:
            #Vdatetimes.append(datetime.datetime.combine(datetime.date(VYearList[ii], VMonthList[ii], VDayList[ii]), j))
            #ii += 1

        MedianSigmaB = gaussian_sigma_to_fwhm * np.array(MedianSigmaB) * pixToArc
        MedianSigmaV = gaussian_sigma_to_fwhm * np.array(MedianSigmaV) * pixToArc
        plt.figure(figsize=(12,8))
        plt.rcParams['font.size'] = 16
        ax = plt.axes()
        ax.scatter(Bdatetimes, MedianSigmaB, label = 'B Filter')
        ax.scatter(Vdatetimes, MedianSigmaV, label = 'V Filter')
        ax.legend()
        ax.set_title("Median FWHM of Stars in Various Exposures")
        ax.set_xlabel('Exposure Start Time', fontsize = 21)
        ax.set_ylabel('Full width half maximum (arcseconds)', fontsize = 23)
        fmt = mdates.DateFormatter('%m-%d %H:%M')
        ax.xaxis.set_major_formatter(fmt)
        #ax.set_xlim(datetime.datetime(2023, 2, 1, 20), datetime.datetime(2023,2,2,3))
        maximum = max(MedianSigmaB) + 0.2
        minimum = min(MedianSigmaV) - 0.2
        ax.set_ylim(minimum,maximum)

        plt.savefig(f'Studentship/Graphs/{len(objectList)}_Objects_ShortSeeing_{dateList[0]}.png')
        plt.show()

    elif mode == 2:
        dateList = dates.split(',')
        objectList = objects.split(',')
        MedianSigmaB2 = []
        MedianSigmaV2 = []
        BYearList = []
        BMonthList = []
        BDayList = []
        VYearList = []
        VMonthList = []
        VDayList = []
        
        
        i = 0
        for i in range(0, len(dateList)):
            fileB = f'Studentship/PhotutilsTutorials/90secExposures{objectList[i]}_B_{dateList[i]}.txt'
            fileV = f'Studentship/PhotutilsTutorials/90secExposures{objectList[i]}_V_{dateList[i]}.txt'
            with open(fileB, 'r') as f1:
                lines = f1.readlines()
            BFile = lines[0].strip('\n')
            with open(fileV, 'r') as f2:
                lines2 = f2.readlines()
            VFile = lines2[0].strip('\n')
            if os.path.exists(f'Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}ReprojectedImages-B.txt') == False:
                imageProject(objectList[i], 'B', BFile, dateList[i], 'single')
                imageProject(objectList[i], 'V', VFile, dateList[i], 'single')
            if method == 'median':
                if (os.path.exists(f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_8_V.fits') == False):
                    imageStack(objectList[i], 'B', 'median', dateList[i], False, 8)
                    imageStack(objectList[i], 'V', 'median', dateList[i], False, 8)
                stackedB = f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_8_B.fits'
                stackedV = f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_8_V.fits'
            else:
                if (os.path.exists(f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_mean_V_8.fits') == False):
                    imageStack(objectList[i], 'B', 'mean', dateList[i], False, 8)
                    imageStack(objectList[i], 'V', 'mean', dateList[i], False, 8)
                stackedB = f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_mean_B_8.fits'
                stackedV = f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_mean_V_8.fits'
            print(objectList[i])
            SigmaListB = SingleImagePSF(stackedB, objectList[i])
            SigmaListV = SingleImagePSF(stackedV, objectList[i])
            MedianSigmaB2.append(np.median(SigmaListB))
            MedianSigmaV2.append(np.median(SigmaListV))
            MedianSigmaB3 = gaussian_sigma_to_fwhm * np.array(MedianSigmaB2) * pixToArc
            MedianSigmaV3 = gaussian_sigma_to_fwhm * np.array(MedianSigmaV2) * pixToArc
            
            headerB = fits.getheader(stackedB, 0)
            startTimeB:str = headerB['DATE-OBS']
            timeStringB = (startTimeB.strip('\'').split('T'))
            timeB = timeStringB[1]
            dateSplitB = (timeStringB[0]).split('-')
            BYearList.append(int(dateSplitB[0]))
            BMonthList.append(int(dateSplitB[1]))
            BDayList.append(int(dateSplitB[2]))   

            headerV = fits.getheader(stackedV, 0)
            startTimeV:str = headerV['DATE-OBS']
            timeStringV = (startTimeV.strip('\'').split('T'))
            timeV = timeStringV[1]
            dateSplitV = (timeStringV[0]).split('-')
            VYearList.append(int(dateSplitV[0]))
            VMonthList.append(int(dateSplitV[1]))
            VDayList.append(int(dateSplitV[2]))   

        Bdates = []
        Vdates = []

        j = 0
        for j in range(0, len(BYearList)):
            Bdates.append(datetime.date(BYearList[j], BMonthList[j], BDayList[j]))
            Vdates.append(datetime.date(VYearList[j], VMonthList[j], VDayList[j]))
        
        plt.figure(figsize=(12,8))
        plt.rcParams['font.size'] = 18
        ax = plt.axes()
        ax.scatter(Bdates, MedianSigmaB3, label = 'B Filter')
        ax.scatter(Vdates, MedianSigmaV3, label = 'V Filter')
        ax.legend()
        ax.set_title("Median FWHM of Stars in Various Exposures")
        ax.set_xlabel('Date', fontsize = 18)
        ax.set_ylabel('Full width half maximum (arcseconds)', fontsize = 18)
        fmt = mdates.DateFormatter('%Y-%b-%d')
        ax.xaxis.set_major_formatter(fmt)
        plt.savefig(f'Studentship/Graphs/{objects}LongSeeing_{dates}.png')
        plt.show()

### SEE IF ITS POSSIBLE TO COLOUR CODE OBJECTS, BUT PUT DIFFERENT SHAPE TOKENS FOR EACH FILTER

'LB1735,WD0830,WD1153,MCT0550,MCT0401'
'20230112,20230112,20230112,20230112,20230112'

plotSeeing('LB1735,WD0830,WD1153,MCT0550,MCT0401','20230112,20230112,20230112,20230112,20230112', 'mean', 1)