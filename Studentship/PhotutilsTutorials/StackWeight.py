import matplotlib.pyplot as plt
from matplotlib.pyplot import*
import numpy as np
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from photutils.background import MMMBackground, MADStdBackgroundRMS
import timeit
from astropy.stats import sigma_clip, sigma_clipped_stats
from ccdproc import combine
from ccdproc import CCDData
from astropy import wcs


def StackWeighting(objectName, filterType, date, stackNum, mode, increment=None, startVal=None):
    """
    increment - only used when mode = 2
    startVal - only used when mode = 2
    the rest of the parameters are controlled by imageStack or ModifiedImageStack
    """

    if mode == 1:
        file = f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{filterType}.txt'
        with open(file, 'r') as f:
            lines = f.readlines()
        
        RMSArray = []
        for line in lines[0:stackNum]:
            image_file = line.strip('\n')
            image_data = fits.getdata(image_file, ext=0)
            image_data = np.array((image_data[~np.isnan(image_data)]))
            bkgrms = MADStdBackgroundRMS()
            image_rms = bkgrms.calc_background_rms(image_data, masked=True)
            RMSArray.append(float(image_rms))

        print(RMSArray)
        
        for i in range(0, len(RMSArray)):
            RMSArray[i] = 1 / RMSArray[i]
        print()
        print(RMSArray)

        sum = np.sum(RMSArray)
        for i in range(0,len(RMSArray)):
            RMSArray[i] = RMSArray[i] / sum
        print()
        print(RMSArray)
        print(np.sum(RMSArray))

    if mode == 2:
        file = f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{filterType}.txt'
        with open(file, 'r') as f:
            lines = f.readlines()
        
        RMSArray = []
        for line in lines[startVal:increment]:
            image_file = line.strip('\n')
            image_data = fits.getdata(image_file, ext=0)
            image_data = image_data[~np.isnan(image_data)]
            bkgrms = MADStdBackgroundRMS()
            image_rms = bkgrms(image_data)
            RMSArray.append(float(image_rms))

        print(RMSArray)
        
        for i in range(0, len(RMSArray)):
            RMSArray[i] = 1 / RMSArray[i]
        print()
        print(RMSArray)

        sum = np.sum(RMSArray)
        for i in range(0,len(RMSArray)):
            RMSArray[i] = RMSArray[i] / sum
        print()
        print(RMSArray)
        print(np.sum(RMSArray))
    return np.array(RMSArray)
    
#StackWeighting('JL163', 'V')

#f3 = open(f'Studentship/PhotutilsTutorials/{objectName}RMSValues-{filterType}.txt', 'a')
#for item in RMSArray:
    #f3.write(str(item) + '\n')
#f3.close()

def weightMean(objectName, FilterType):
    file = 'Studentship/PhotutilsTutorials/JL163ReprojectedImages-V.txt'
    with open(file, 'r') as f:
        lines = f.readlines()
        
    CCDObjects = []
    for line in lines:
        image_file = line.strip('\n')
        image_data = fits.getdata(image_file, ext=0)
        hdu1 = fits.open(image_file)[0]     
        w = wcs.WCS(hdu1.header)
        CCDObjects.append(CCDData(image_data, wcs = w, unit = 'adu'))

    weighted_median = combine(CCDObjects, method = 'average', weights = StackWeighting('JL163', 'V'), sigma_clip = True, sigma_clip_high_thresh=3.0, sigma_clip_low_thresh=3.0)
    data = weighted_median.data
    weighted_median.header = hdu1.header
    norm = ImageNormalize(stretch = LogStretch())
    print(weighted_median.header)

    #fits.writeto(f'Studentship/ReprojectedImages/{objectName}Reprojected-{num}-{filterType}.fits', data, weighted_median.header, overwrite=True)
    plt.imshow(data, 'gray', norm = norm)
    plt.show()



