import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from photutils.background import MMMBackground, MADStdBackgroundRMS
from scipy.stats import linregress

def NoiseGraph(objectName, date:str, method, stackNum, increment, file1, file2):
    """
    
    date in format - yyyymmdd,yyyymmdd,yyyymmdd.... for however many datasets were used
    method - 'mean' or 'median'
    objectName - self explanatory, often shortened and with no hyphens. e.g JL163 instead of MCT-0401-4017
    stackNum - total number of images used in largest stack
    increment - the difference in number of images between successive stacked images
    file1 and file2 - same as in FluxCalibrationAnalysis. File1 is the original B image, File2 is the original V image
    
    """
    originalIncrement = increment
    BFiles = [file1]
    VFiles = [file2]
    time = [90]
    exposure = 90
    zeroPointBList = []
    zeroPointVList = []


    #THIS CODE IS REPEATED TO GET ZERO POINT CONSTANTS FOR THE FIRST IMAGES, WHICH HAVE NOT BEEN INDEPENDENTLY ANALYSED
    file1 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_{increment}_{method}-Coefficients.txt'
    with open(file1, 'r') as f:
        lines = f.readlines()
        var1 = lines[3].strip('\n').split(', ')
        var2 = lines[4].strip('\n').split(', ') 
        zeroPointBList.append(float(var1[1]))
        zeroPointVList.append(float(var2[1]))

    while increment <= stackNum:
        file1 = f'Studentship/PhotutilsTutorials/{date}_{objectName}_{increment}_{method}-Coefficients.txt'
        with open(file1, 'r') as f:
            lines = f.readlines()
            var1 = lines[3].strip('\n').split(', ')
            var2 = lines[4].strip('\n').split(', ')
            zeroPointBList.append(float(var1[1]))
            zeroPointVList.append(float(var2[1]))
        if method == 'median':
            BFiles.append(f'Studentship/StackedImages/{date}_{objectName}Stacked_{increment}_B.fits')
            VFiles.append(f'Studentship/StackedImages/{date}_{objectName}Stacked_{increment}_V.fits')
        elif method == 'mean':
            BFiles.append(f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_B_{increment}.fits')
            VFiles.append(f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_V_{increment}.fits')
        newTime = exposure * increment
        time.append(newTime)
        increment += originalIncrement
    
    i = 0
    BNoiseList = []
    VNoiseList = []
    bkgrms = MADStdBackgroundRMS()
    for i in range(0, len(BFiles)):
        B_image = fits.getdata(BFiles[i], ext=0)
        V_image = fits.getdata(VFiles[i], ext=0)
        B_image = B_image[100:700, 480:1230]
        V_image = V_image[100:700, 480:1230]
        B_rms = bkgrms(B_image)
        V_rms = bkgrms(V_image)
        BNoiseList.append(B_rms)
        VNoiseList.append(V_rms)


    instrumentalMagListB = []
    instrumentalMagListV = []
    jj = 0
    for jj in range(0, len(BNoiseList)):
        instrumentalMagListB.append(-2.5*np.log10(float(BNoiseList[jj] / 90)))
        instrumentalMagListV.append(-2.5*np.log10(float(VNoiseList[jj] / 90)))
    ii = 0
    for ii in range(0, len(instrumentalMagListB)):
        instrumentalMagListB[ii] = instrumentalMagListB[ii] + zeroPointBList[ii]
        instrumentalMagListV[ii] = instrumentalMagListV[ii] + zeroPointVList[ii]
    
    magArcsecondListB = []
    magArcsecondListV = []
    ii = 0
    for ii in range(0, len(instrumentalMagListB)):
        magArcsecondListB.append(instrumentalMagListB[ii] + 2.5*np.log10(0.8305**2))
        magArcsecondListV.append(instrumentalMagListV[ii] + 2.5*np.log10(0.8305**2))

    logTime = []
    for item in time:
        logTime.append(np.log(item))
    regressB = linregress(logTime[3:], magArcsecondListB[3:])
    regressV = linregress(logTime[3:], magArcsecondListV[3:])
    slopeB = regressB.slope
    interceptB = regressB.intercept
    print(slopeB)

    slopeV = regressV.slope
    interceptV = regressV.intercept
    print(slopeV)
    y1 = magArcsecondListB[1]
    y2 = magArcsecondListB[2]
    x1 = time[1]
    x2 = time[2]
    a1 = float((y1-y2) / (np.sqrt(x1) - np.sqrt(x2)))
    b1 = y2 - a1*np.sqrt(x2)
    x = np.linspace(0,3000)
    y = a1 * np.sqrt(x) + b1


    y3 = magArcsecondListV[1]
    y4 = magArcsecondListV[2]
    x3 = time[1]
    x4 = time[2]
    a2 = float((y3-y4) / (np.sqrt(x3) - np.sqrt(x4)))
    b2 = y4 - a2*np.sqrt(x4)
    y5 = a2 * np.sqrt(x) + b2
    
    plt.figure(figsize=(12,8))
    plt.rcParams['font.size'] = 16
    plt.scatter(time, magArcsecondListB, label = 'B Filter')
    plt.scatter(time, magArcsecondListV, label = 'V Filter')
    plt.plot(x,y, label = 'Theoretical Square Root, B Filter')
    plt.plot(x,y5, label = 'Theoretical Square Root, V Filter')
    plt.legend()
    plt.title(f'Noise vs Time, {objectName}')
    plt.xlabel('Time (seconds)')
    plt.ylabel("Noise (magnitudes per arcsecond squared)")
    plt.xlim(0,5250)
    plt.ylim(25.5, 22.5)
    plt.savefig(f'Studentship/Graphs/{objectName}_SubtractedNoiseVTime_{method}.png')
    plt.show()

    x = np.linspace(6.5,np.log(5000))
    y = slopeB * x + interceptB
    y2 = slopeV * x + interceptV
    plt.scatter(logTime, magArcsecondListB, label = 'B Filter')
    plt.scatter(logTime, magArcsecondListV, label = 'V Filter')
    plt.legend()
    plt.plot(x,y)
    plt.plot(x,y2)
    plt.title(f'Noise vs ln(Time), {objectName}')
    plt.xlabel('ln(Time)')
    plt.ylabel("Noise (pixels)")
    plt.show()
    
    

NoiseGraph('MCT0401', '20230112,20230114,20230118,20230126,20230127', 'median', 40, 4, 'Studentship/SPIRIT6-Observations/GoodDarks/Multiple-Exposures/20230112/MCT-0401-4017/MCT-0401-4017-S001-B-R001-B.fts', 'Studentship/SPIRIT6-Observations/GoodDarks/Multiple-Exposures/20230112/MCT-0401-4017/MCT-0401-4017-S001-V-R001-V.fts')