
from distutils.log import Log
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clip, sigma_clipped_stats
from ccdproc import combine
from ccdproc import CCDData
from astropy import wcs
from StackWeight import StackWeighting
from photutils.background import MMMBackground, MADStdBackgroundRMS
#objectName = 'JL163'
#FilterType = 'V'

def imageStack(objectName, FilterType, method, date, bkgsub, stackNum=None):
    """
    objectName - the name of an individual object, e.g. 'JL163' \n
    filterType - letter corresponding to filter, 'B' or 'V' \n
    originalFileName - this is the base image that all other images will be projected onto. give the path to the image. \n
    method - 'median' or 'mean' \n
    date - yyyymmdd format e.g. '20230112' for Jan 12th 2023. \n
    stackNum - how many images from the total dataset you want to stack together.
    \n 
    
    IF mode from FluxCalibrationAnalysis is multiple: \n
        date - yyyymmdd,yyyymmdd.... repeated for as many datasets as you want of the given object e.g. 20230112,20230114,20230118... \n
    
    """
    fileNames = []
    file = f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{FilterType}.txt'
    with open(file, 'r') as f1:
        lines = f1.readlines()
    i = 0
    for line in lines:
        i += 1
    if stackNum == None:
        stackNum = i

    if method == 'median':
        with open(file, 'r') as f:
            lines = f.readlines()
            num = np.shape(lines)[0]
            initFile = lines[0].strip('\n')
            initHeader = fits.open(initFile)[0]
            initData = fits.getdata(initFile, ext=0)
            i = 1
            dataList = []
            hduList = []
            dataList.append(initData)
            hduList.append(initHeader)
            for line in lines[1:stackNum]:
                line1 = line.strip('\n')
                hdu = fits.open(line1)[0]
                data = fits.getdata(line1, ext=0)
                if bkgsub == True:
                    mmm_bkg = MMMBackground()
                    img_bkg = mmm_bkg(data)
                    bkg_array = img_bkg * np.ones(shape=np.shape(data))
                    data = np.subtract(data,bkg_array)
                dataList.append(data)
                hduList.append(hdu)
            d = np.stack(dataList)
            final = np.median(d, axis=0)
            fits.writeto(f'Studentship/StackedImages/{date}_{objectName}Stacked_{stackNum}_{FilterType}.fits', final, initHeader.header, overwrite=True)
            fileNames.append(f'Studentship/StackedImages/{date}_{objectName}Stacked_{stackNum}_{FilterType}.fits')
            

    
        f2 = open(f'Studentship/PhotutilsTutorials/{date}_{objectName}MedianStackedImages-{FilterType}.txt', 'a')
        f2.write(initFile + '\n')
        for item in fileNames:
            f2.write(item + '\n')
        f2.close()
    elif method == 'mean':
        file4 = f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{FilterType}.txt'
        with open(file4, 'r') as f4:
            lines4 = f4.readlines()
            
        CCDObjects = []
        jj = 1
        for line in lines4:
            if jj <= stackNum:
                image_file = line.strip('\n')
                image_data = fits.getdata(image_file, ext=0)
                if bkgsub == True:
                    mmm_bkg = MMMBackground()
                    img_bkg = mmm_bkg(image_data)
                    bkg_array = img_bkg * np.ones(shape=np.shape(image_data))
                    image_data = np.subtract(image_data, bkg_array)
                hdu1 = fits.open(image_file)[0]     
                w = wcs.WCS(hdu1.header)
                CCDObjects.append(CCDData(image_data, wcs = w, unit = 'adu'))
            jj += 1
        print(len(CCDObjects))

        weighted_mean = combine(CCDObjects, method = 'average', weights = StackWeighting(objectName, FilterType, date, stackNum, 1), sigma_clip = True, sigma_clip_high_thresh=2.0, sigma_clip_low_thresh=2.0)
        data = weighted_mean.data
        weighted_mean.header = hdu1.header
        fits.writeto(f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_{FilterType}_{stackNum}.fits', data, weighted_mean.header, overwrite=True)
        f2 = open(f'Studentship/PhotutilsTutorials/{date}_{objectName}MeanStackedImages-{FilterType}_{stackNum}.txt', 'a')
        f2.write( f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_{FilterType}_{stackNum}.fits' + '\n')
        f2.close()




def ModifiedImageStack(objectName, FilterType, method, date, stackNum=None, increment=None):
    """
    objectName - the name of an individual object, e.g. 'JL163' \n
    filterType - letter corresponding to filter, 'B' or 'V' \n
    originalFileName - this is the base image that all other images will be projected onto. give the path to the image. \n
    method - 'median' or 'mean' \n
    date - yyyymmdd format e.g. '20230112' for Jan 12th 2023. \n
    stackNum - how many images from the total dataset you want to stack together.\n
    increment - almost always will be 4. The depth of stacking for each stacked image.\n
    \n 
    
    IF mode from FluxCalibrationAnalysis is multiple: \n
        date - yyyymmdd,yyyymmdd.... repeated for as many datasets as you want of the given object e.g. 20230112,20230114,20230118... \n
    
    """
    fileNames = []
    file = f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{FilterType}.txt'
    with open(file, 'r') as f1:
        lines = f1.readlines()
    i = 0
    for line in lines:
        i += 1
    if stackNum == None:
        stackNum = i

    if method == 'median':
        with open(file, 'r') as f:
            lines = f.readlines()
            num = np.shape(lines)[0]
            initFile = lines[0].strip('\n')
            initHeader = fits.open(initFile)[0]
            initData = fits.getdata(initFile, ext=0)
            originalIncrement = increment
            startVal = 0
            while increment <= stackNum:
                dataList = []
                hduList = []
                for line in lines[startVal:increment]:
                    line1 = line.strip('\n')
                    hdu = fits.open(line1)[0]
                    data = fits.getdata(line1, ext=0)
                    dataList.append(data)
                    hduList.append(hdu)
                d = np.stack(dataList)
                final = np.median(d, axis=0)
                fits.writeto(f'Studentship/StackedImages/{date}_{objectName}Stacked_{startVal+1}-{increment}_{FilterType}.fits', final, initHeader.header, overwrite=True)
                fileNames.append(f'Studentship/StackedImages/{date}_{objectName}Stacked_{startVal+1}-{increment}_{FilterType}.fits')
                increment += originalIncrement
                startVal += originalIncrement

    elif method == 'mean':
        file4 = f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{FilterType}.txt'
        with open(file4, 'r') as f4:
            lines4 = f4.readlines()
        originalIncrement = increment
        startVal = 0
        while increment <= stackNum:
            CCDObjects = []
            jj = 1
            for line in lines4[startVal:increment]:
                if jj <= stackNum:
                    image_file = line.strip('\n')
                    image_data = fits.getdata(image_file, ext=0)
                    hdu1 = fits.open(image_file)[0]     
                    w = wcs.WCS(hdu1.header)
                    CCDObjects.append(CCDData(image_data, wcs = w, unit = 'adu'))
                jj += 1

            weighted_mean = combine(CCDObjects, method = 'average', weights = StackWeighting(objectName, FilterType, date, stackNum, 2, increment, startVal), sigma_clip = True, sigma_clip_high_thresh=2.0, sigma_clip_low_thresh=2.0)
            data = weighted_mean.data
            weighted_mean.header = hdu1.header
            fits.writeto(f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_{FilterType}_{startVal+1}-{increment}.fits', data, weighted_mean.header, overwrite=True)
            startVal += originalIncrement
            increment += originalIncrement


#imageStack('WD1153', 'V', 'mean', '20230112,20230114,20230118', 24)


