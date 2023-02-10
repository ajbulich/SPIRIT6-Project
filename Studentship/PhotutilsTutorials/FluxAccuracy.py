
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.io import fits
from FluxCalibration1 import FluxCalibration
from projectMultiple import imageProject
from stackMultiple import imageStack
from FinalAperturePhotometry import* 
from psfForStackedImages import*
from astropy.stats import sigma_clip, sigma_clipped_stats, gaussian_sigma_to_fwhm
import os
from psfForStackedImages import SingleImagePSF
import datetime
import matplotlib.dates as mdates
import timeit

def FluxAccuracy(objects:str, dates:str, method):

    """
    objects - string containing all objects, including repetition for different dates. e.g. 'JL163,JL163,WD0830,WD0830,WD1153...' \n
    dates - corresponding dates for those objects. e.g. 20230127,20230125,20230118... Must be same length as objects \n
    method - median or mean
    """
    objectList = objects.split(',')
    dateList = dates.split(',')
    i = 0
    diffListB = []
    diffListV = []
    YearList = []
    MonthList = []
    DayList = []
    for i in range(0, len(dateList)):
        fileB = f'Studentship/PhotutilsTutorials/90secExposures{objectList[i]}_B_{dateList[i]}.txt'
        fileV = f'Studentship/PhotutilsTutorials/90secExposures{objectList[i]}_V_{dateList[i]}.txt'
        with open(fileB, 'r') as f1:
            lines = f1.readlines()
        BFile = lines[0].strip('\n')
        with open(fileV, 'r') as f2:
            lines2 = f2.readlines()
        VFile = lines2[0].strip('\n')
        if (os.path.exists(f'Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}ReprojectedImages-B.txt') == False) or (os.path.exists(f'Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}ReprojectedImages-V.txt') == False):
            if os.path.exists(f'Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}ReprojectedImages-B.txt'):
                os.remove(f'Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}ReprojectedImages-B.txt')
            if os.path.exists(f'Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}ReprojectedImages-V.txt'):
                os.remove(f'Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}ReprojectedImages-V.txt')
            imageProject(objectList[i], 'B', BFile, dateList[i], 'single')
            imageProject(objectList[i], 'V', VFile, dateList[i], 'single')
        if method == 'median':
            if (os.path.exists(f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_8_V.fits') == False):
                imageStack(objectList[i], 'B', 'median', dateList[i], False, 8)
                imageStack(objectList[i], 'V', 'median', dateList[i], False, 8)
            if os.path.exists(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}Stacked_B_8_PSF.txt"):
                os.remove(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}Stacked_B_8_PSF.txt")
                os.remove(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}Stacked_V_8_PSF.txt")
                os.remove(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}_MedianStackedApertureResults-B_8.txt")
                os.remove(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}_MedianStackedApertureResults-V_8.txt")
            stackedB = f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_8_B.fits'
            stackedV = f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_8_V.fits'
    
        else:
            if (os.path.exists(f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_mean_V_8.fits') == False):
                imageStack(objectList[i], 'B', 'mean', dateList[i], False, 8)
                imageStack(objectList[i], 'V', 'mean', dateList[i], False, 8)
            if (os.path.exists(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}Stacked_B_8_mean_PSF.txt")):
                os.remove(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}Stacked_B_8_mean_PSF.txt")
                os.remove(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}Stacked_V_8_mean_PSF.txt")
                os.remove(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}_MeanStackedApertureResults-B_8.txt")
                os.remove(f"Studentship/PhotutilsTutorials/{dateList[i]}_{objectList[i]}_MeanStackedApertureResults-V_8.txt")
            stackedB = f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_mean_B_8.fits'
            stackedV = f'Studentship/StackedImages/{dateList[i]}_{objectList[i]}Stacked_mean_V_8.fits'
            
        PSFStackedPhotometry(1, 'B', objectList[i], method, dateList[i], 8)
        PSFStackedPhotometry(1, 'V', objectList[i], method, dateList[i], 8)
        annulusPhotometry(objectList[i], 'B', method, dateList[i], 8)
        annulusPhotometry(objectList[i], 'V', method, dateList[i], 8)
        Diff = FluxCalibration(objectList[i], method, dateList[i], 8)
        BDiff = Diff[0]
        VDiff = Diff[1]
        diffListB.append(np.abs(np.median(BDiff)))
        diffListV.append(np.abs(np.median(VDiff)))

        headerB = fits.getheader(stackedB, 0)
        startTimeB:str = headerB['DATE-OBS']
        timeStringB = (startTimeB.strip('\'').split('T'))
        timeB = timeStringB[1]
        dateSplitB = (timeStringB[0]).split('-')
        YearList.append(int(dateSplitB[0]))
        MonthList.append(int(dateSplitB[1]))
        DayList.append(int(dateSplitB[2]))   
    if method == 'mean':
        f3 = open('Studentship/MagnitudeResults/MeanCollectiveMagnitudeResultsB.txt', 'a')
        i = 0
        for i in range(0, len(diffListB)):
            f3.write(str(diffListB[i]) + '\n')
            i += 1
        f3.close()

        f4 = open('Studentship/MagnitudeResults/MeanCollectiveMagnitudeResultsV.txt', 'a')
        i = 0
        for i in range(0, len(diffListB)):
            f4.write(str(diffListV[i]) + '\n')
            i += 1
        f4.close()

    elif method == 'median':
        f3 = open('Studentship/MagnitudeResults/MedianCollectiveMagnitudeResultsB.txt', 'a')
        i = 0
        for i in range(0, len(diffListB)):
            f3.write(str(diffListB[i]) + '\n')
            i += 1
        f3.close()

        f4 = open('Studentship/MagnitudeResults/MedianCollectiveMagnitudeResultsV.txt', 'a')
        i = 0
        for i in range(0, len(diffListB)):
            f4.write(str(diffListV[i]) + '\n')
            i += 1
        f4.close()


    j = 0
    Dates = []
    for j in range(0, len(YearList)):
        Dates.append(datetime.date(YearList[j], MonthList[j], DayList[j]))

    plt.figure(figsize=(12,8))
    ax = plt.axes()
    ax.scatter(Dates, diffListB, label = 'B Filter')
    ax.scatter(Dates, diffListV, label = 'V Filter')
    ax.legend()
    ax.set_title("Flux Calibration Accuracy")
    ax.set_xlabel('Date')
    ax.set_ylabel('Median Magnitude Error (magnitude/pixel)')
    fmt = mdates.DateFormatter('%Y-%b-%d')
    ax.xaxis.set_major_formatter(fmt)
    plt.savefig(f'Studentship/Graphs/{objects}_{method}_MagError.png')
    plt.show()


startTime = timeit.default_timer()
FluxAccuracy('JL163,JL163,JL163,WD0830,WD0830,WD0830,WD0830,WD0830,WD0830,WD0830,WD1153,WD1153,WD1153,WD1153,WD1153,WD1153', '20221214,20230115,20230116,20221214,20230112,20230114,20230118,20230125,20230127,20230128,20230112,20230114,20230118,20230125,20230127,20230128', 'mean')
endTime = timeit.default_timer()
runTime = endTime - startTime
print("The total runtime was: ", runTime)


