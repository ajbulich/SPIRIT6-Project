import matplotlib.pyplot as plt
from matplotlib.pyplot import*
import numpy as np
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from photutils.background import MMMBackground, MADStdBackgroundRMS
from stackMultiple import*

def HistogramAnalysis(objectName, date, FilterType, stackNum, increment, method, mode):

    """
    objectName - single string for object name \n
    date - string containing date or dates in yyyymmdd,yyyymmdd.... format \n
    FilterType - 'B' or 'V'\n
    stackNum - maximum size of stack\n
    increment - increment that you wish to create and analyse the stacked images in. \n
                if increment is 4, every 4 images a new stack is created and a new histogram made\n
    method - 'mean' or 'median' \n
    MODE - 1 is for normal histograms with cumulative images\n
    MODE - 2 is for every 4 images, create a stack and then histogram those\n
    """
    if mode == 1:
        originalIncrement = increment
        while increment <= stackNum:
            if method == 'median':
                image_file = f'Studentship/StackedImages/{date}_{objectName}Stacked_{increment}_{FilterType}.fits'
            elif method == 'mean':
                image_file = f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_{FilterType}_{increment}.fits'

            image_data = fits.getdata(image_file, ext=0)[400:1200, 680:1380]
            counts3, bins3 = np.histogram(image_data[~np.isnan(image_data)], bins=32000)
            mmm_bkg = MMMBackground()
            bkg1 = mmm_bkg(image_data[~np.isnan(image_data)])
            bkgrms = MADStdBackgroundRMS()
            rms1 = bkgrms(image_data[~np.isnan(image_data)])
            print(image_file)
            plt.figure(figsize = (16,8))
            plt.hist(bins3[:-1], bins3[0:200], weights = counts3, alpha = 0.5)
            plt.axvline(bkg1, color = 'blue')
            plt.axvline(bkg1 - rms1, color = 'blue', ls = '--', alpha = 0.5)
            plt.axvline(bkg1 + rms1, color = 'blue', ls = '--', alpha = 0.5)
            plt.yscale('log')
            plt.title("Histogram")
            plt.savefig(f'Studentship/Graphs/{objectName}_{date}_{method}_{FilterType}_SubtractedNoiseHistogram_{increment}.png')
            plt.show()
            

            increment += originalIncrement

    if mode == 2:
        ModifiedImageStack(objectName, FilterType, method, date, stackNum, increment)
        originalIncrement = increment
        startVal = 0
        while increment <= stackNum:
            if method == 'median':
                image_file = f'Studentship/StackedImages/{date}_{objectName}Stacked_{startVal+1}-{increment}_{FilterType}.fits'
            elif method == 'mean':
                image_file = f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_{FilterType}_{startVal+1}-{increment}.fits'

            image_data = fits.getdata(image_file, ext=0)[400:1200, 680:1380]
            counts3, bins3 = np.histogram(image_data[~np.isnan(image_data)], bins=32000)
            mmm_bkg = MMMBackground()
            bkg1 = mmm_bkg(image_data[~np.isnan(image_data)])
            bkgrms = MADStdBackgroundRMS()
            rms1 = bkgrms(image_data[~np.isnan(image_data)])
            print(image_file)
            plt.figure(figsize = (16,8))
            plt.hist(bins3[:-1], bins3[0:120], weights = counts3, alpha = 0.5)
            plt.axvline(bkg1, color = 'blue')
            plt.axvline(bkg1 - rms1, color = 'blue', ls = '--', alpha = 0.5)
            plt.axvline(bkg1 + rms1, color = 'blue', ls = '--', alpha = 0.5)
            plt.yscale('log')
            plt.title("Histogram")
            plt.savefig(f'Studentship/Graphs/{objectName}_{date}_{method}_{FilterType}_SubtractedNoiseHistogram_{startVal}-{increment}.png')
            plt.show()
            

            increment += originalIncrement
            startVal += originalIncrement

HistogramAnalysis("WD0830", '20221214,20230112,20230114,20230118', 'V', 32, 4, 'mean', 1)