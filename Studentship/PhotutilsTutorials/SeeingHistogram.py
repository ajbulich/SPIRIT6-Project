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
pixToArc = 0.8305

def addData(mode, dates:str, objects:str=None):   
    """
    mode = 1\n
        - dates will be a single date, objects will be a string of csv's representing all objects imaged that night\n
          adds all data to a file for that specific night.\n
    \n
    mode = 2\n
        - unique setup. Have to run mode 1 for every single date first. Then subsequently we run mode 2.\n
          The idea is that it will access all the text files for every date, and append all the data to one\n
          dataset and then plot this.\n
          dates will be yyyymmdd,yyyymmdd,yyyymmdd,..... for all dates we've used mode1 for so far.\n
    """
    if mode == 1:
        dataListB = []
        dataListV = []
        objectList = objects.split(',')
        for item in objectList:
            file = f'Studentship/PhotutilsTutorials/90secExposures{item}_B_{dates}.txt'
            file2 = f'Studentship/PhotutilsTutorials/90secExposures{item}_V_{dates}.txt'
            with open(file, 'r') as f:
                lines = f.readlines()
            for line in lines:
                line = line.strip('\n')
                data = gaussian_sigma_to_fwhm * pixToArc * np.array(SingleImagePSF(line, item))
                for item2 in data:
                    dataListB.append(item2)
            with open(file2, 'r') as f2:
                lines2 = f2.readlines()
            for line2 in lines2:
                line2 = line2.strip('\n')
                data2 = gaussian_sigma_to_fwhm * pixToArc * np.array(SingleImagePSF(line2, item))
                for item3 in data2:
                    dataListV.append(item3)
        
        f3 = open(f'Studentship/PhotutilsTutorials/FWHMList_{dates}_B.txt', 'a')
        for item4 in dataListB:
            f3.write(str(item4)+'\n')
        f3.close()

        f4 = open(f'Studentship/PhotutilsTutorials/FWHMList_{dates}_V.txt', 'a')
        for item5 in dataListV:
            f4.write(str(item5)+'\n')
        f4.close()

    elif mode == 2:
        dateList = dates.split(',')
        dataListB2 = []
        dataListV2 = []
        for item in dateList:
            Bfile = f'Studentship/PhotutilsTutorials/FWHMList_{item}_B.txt'
            Vfile = f'Studentship/PhotutilsTutorials/FWHMList_{item}_V.txt'

            with open(Bfile, 'r') as f5:
                lines5 = f5.readlines()
                for line5 in lines5:
                    line5 = float(line5.strip('\n'))
                    dataListB2.append(line5)

            with open(Vfile, 'r') as f6:
                lines6 = f6.readlines()
                for line6 in lines6:
                    line6 = float(line6.strip('\n'))
                    dataListV2.append(line6)

        f7 = open('Studentship/PhotutilsTutorials/AllFWHM_B.txt', 'a')
        for item in dataListB2:
            f7.write(str(item)+'\n')
        f7.close()

        f8 = open('Studentship/PhotutilsTutorials/AllFWHM_V.txt', 'a')
        for item2 in dataListV2:
            f8.write(str(item2)+'\n')
        f8.close()


def createHistogram(mode, date=None):
    """
    if mode = 1:\n
        must specify a date\n
    elif mode = 2:\n
        don't have to specify a date\n
    """
    if mode == 1:  
        dataListB = []
        dataListV = []
        Bfile = f'Studentship/PhotutilsTutorials/FWHMList_{date}_B.txt'
        Vfile = f'Studentship/PhotutilsTutorials/FWHMList_{date}_V.txt'
        with open(Bfile, 'r') as f:
            lines = f.readlines()
        for line in lines:
            line = float(line.strip('\n'))
            dataListB.append(line)
        with open(Vfile, 'r') as f2:
            lines2 = f2.readlines()
        for line2 in lines2:
            line2 = float(line2.strip('\n'))
            dataListV.append(line2)
        
        dataListV = sigma_clip(dataListV, sigma=6, cenfunc='median', stdfunc = 'mad_std', axis=0, masked=False)
        dataListV = dataListV[~np.isnan(dataListV)]
        dataListB = sigma_clip(dataListB, sigma=6, cenfunc='median', stdfunc = 'mad_std', axis=0, masked=False)
        dataListB = dataListB[~np.isnan(dataListB)]
        counts, bins = np.histogram(dataListB, bins=24)
        counts2, bins2 = np.histogram(dataListV, bins=24)
        bkg = MMMBackground()
        rms = MADStdBackgroundRMS()
        B_bkg = bkg(dataListB)
        B_rms = rms(dataListB)
        V_bkg = bkg(dataListV)
        V_rms = rms(dataListV)
        plt.figure(figsize = (16,8))
        plt.rcParams['font.size'] = 16
        plt.hist(bins[:-1], bins, weights = counts, alpha = 0.5, label = 'B Filter')
        plt.hist(bins2[:-1], bins2, weights = counts2, alpha = 0.5, label = 'V Filter')
        plt.axvline(B_bkg, color = 'blue',  label = 'MMM B')
        plt.axvline(B_bkg - B_rms, color = 'blue', ls = '--', alpha = 0.5)
        plt.axvline(B_bkg + B_rms, color = 'blue', ls = '--', alpha = 0.5)
        plt.axvline(V_bkg, color = 'red', label = 'MMM V')
        plt.axvline(V_bkg - V_rms, color = 'red', ls = '--', alpha = 0.5)
        plt.axvline(V_bkg + V_rms, color = 'red', ls = '--', alpha = 0.5)
        plt.xlim(0)
        plt.xlabel("FWHM Values (arcseconds)", fontsize = 21)
        plt.title(f'{date} FWHM Histogram')
        plt.legend()
        plt.savefig(f"Studentship/Graphs/FWHM_Histogram_{date}.png")
        plt.show()

    elif mode == 2:
        Bfile = 'Studentship/PhotutilsTutorials/AllFWHM_B.txt'
        Vfile = 'Studentship/PhotutilsTutorials/AllFWHM_V.txt'
        dataListB2 = []
        dataListV2 = []
        with open(Bfile, 'r') as f3:
            lines3 = f3.readlines()
        for line3 in lines3:
            line3 = float(line3.strip('\n'))
            dataListB2.append(line3)
        with open(Vfile, 'r') as f4:
            lines4 = f4.readlines()
        for line4 in lines4:
            line4 = float(line4.strip('\n'))
            dataListV2.append(line4)

        dataListV2 = sigma_clip(dataListV2, sigma=8, cenfunc='median', stdfunc = 'mad_std', axis=0, masked=False)
        dataListV2 = dataListV2[~np.isnan(dataListV2)]
        dataListB2 = sigma_clip(dataListB2, sigma=8, cenfunc='median', stdfunc = 'mad_std', axis=0, masked=False)
        dataListB2 = dataListB2[~np.isnan(dataListB2)]

        counts, bins = np.histogram(dataListB2, bins=32)
        counts2, bins2 = np.histogram(dataListV2, bins=32)
        bkg = MMMBackground()
        rms = MADStdBackgroundRMS()
        B_bkg = bkg(dataListB2)
        B_rms = rms(dataListB2)
        V_bkg = bkg(dataListV2)
        V_rms = rms(dataListV2)
        plt.figure(figsize = (16,8))
        plt.rcParams['font.size'] = 14
        plt.hist(bins[:-1], bins[1:], weights = counts, alpha = 0.5, label = 'B Filter')
        plt.hist(bins2[:-1], bins2[1:], weights = counts2, alpha = 0.5, label = 'V Filter')
        plt.axvline(B_bkg, color = 'blue', label='MMMBackground B Filter')
        plt.axvline(B_bkg - B_rms, color = 'blue', ls = '--', alpha = 0.5, label='MAD B Filter')
        plt.axvline(B_bkg + B_rms, color = 'blue', ls = '--', alpha = 0.5)
        plt.axvline(V_bkg, color = 'red', label = 'MMMBackground V Filter')
        plt.axvline(V_bkg - V_rms, color = 'red', ls = '--', alpha = 0.5, label='MAD V Filter')
        plt.axvline(V_bkg + V_rms, color = 'red', ls = '--', alpha = 0.5)
        plt.xlabel("FWHM Values (arcseconds)", fontsize = 18)
        plt.legend()
        plt.xlim(0)
        plt.title("Seeing in Perth, all data")
        plt.savefig("Studentship/Graphs/AllFWHM_Histogram.png")
        plt.show()
        print(f'{B_bkg} +- {B_rms}')
        print(f'{V_bkg} +- {V_rms}')

#addData(1, '20230202-(3)', 'LB1735,MCT0550,WD0830,WD1153') #DO NOT INCLUDE 20230131 in the final histogram. Its munted.
createHistogram(1, '20230112')      


        

