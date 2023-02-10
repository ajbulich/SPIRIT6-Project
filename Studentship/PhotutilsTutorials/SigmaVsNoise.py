import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.io import fits
from stackMultiple import imageStack
from FinalAperturePhotometry import* 
from psfForStackedImages import*
from astropy.stats import sigma_clip, sigma_clipped_stats
import os

def SigmaVsNoise(objectName, date, method, FilterType, stackNum, increment):
    """
    makes a graph of sigma vs noise
    """
    noiseList = []
    sigmaList = []
    originalIncrement = increment
    while increment <= stackNum:
        if method == 'median':
            image_file = f'Studentship/StackedImages/{date}_{objectName}Stacked_{increment}_{FilterType}.fits'
            PSF_file = f"Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_{FilterType}_{increment}_PSF.txt"
        elif method == 'mean':
            image_file = f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_{FilterType}_{increment}.fits'
            PSF_file = f"Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_{FilterType}_{increment}_mean_PSF.txt"
        rms = MADStdBackgroundRMS()
        image_data = fits.getdata(image_file, ext=0)
        image_rms = rms(image_data)
        noiseList.append(image_rms)
        with open(PSF_file, 'r') as f:
            lines = f.readlines()
        sigmaList2 = []
        for line in lines:
            var = line.split(' ')
            sigmaList2.append(float(var[3]))
        medianSigma = np.median(sigmaList2)
        sigmaList.append(medianSigma)

        increment += originalIncrement

    plt.scatter(sigmaList, noiseList)
    plt.xlabel('sigma (pixels)')
    plt.ylabel('noise (pixels)')
    plt.show()

    
        
SigmaVsNoise('WD0830', '20221214,20230112,20230118,20230114', 'median', 'V', 32, 4)


        
