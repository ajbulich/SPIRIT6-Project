import numpy as np
import matplotlib.pyplot as plt
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.stats import sigma_clip


def createFluxHistogram(method):
    dataListB = []
    dataListV = []
    if method == 'median':
        fileB = 'Studentship/MagnitudeResults/MedianCollectiveMagnitudeResultsB.txt'
        fileV = 'Studentship/MagnitudeResults/MedianCollectiveMagnitudeResultsV.txt'
    elif method == 'mean':
        fileB = 'Studentship/MagnitudeResults/MeanCollectiveMagnitudeResultsB.txt'
        fileV = 'Studentship/MagnitudeResults/MeanCollectiveMagnitudeResultsV.txt'
    with open(fileB, 'r') as f:
        linesB = f.readlines()

    for line in linesB:
        line = float(line.strip('\n'))
        dataListB.append(line)

    with open(fileV, 'r') as f2:
        linesV = f2.readlines()

    for line2 in linesV:
        line2 = float(line2.strip('\n'))
        dataListV.append(line2)

    #dataListB = sigma_clip(dataListB, sigma=10, cenfunc='median', stdfunc = 'mad_std', axis=0, masked=False)
    #dataListB = dataListB[~np.isnan(dataListB)]
    dataListV = sigma_clip(dataListV, sigma=3, cenfunc='median', stdfunc = 'mad_std', axis=0, masked=False)
    dataListV = dataListV[~np.isnan(dataListV)]
    counts, bins = np.histogram(dataListB, bins=24)
    counts2, bins2 = np.histogram(dataListV, bins=16)

    bkg = MMMBackground()
    rms = MADStdBackgroundRMS()
    B_bkg = bkg(dataListB)
    B_rms = rms(dataListB)
    V_bkg = bkg(dataListV)
    V_rms = rms(dataListV)
    plt.figure(figsize = (16,8))
    plt.rcParams['font.size'] = 14
    plt.hist(bins[:-1], bins, weights = counts, alpha = 0.5, label = 'B Filter')
    plt.hist(bins2[:-1], bins2, weights = counts2, alpha = 0.5, label = 'V Filter')
    plt.xlabel("Absolute Error in Star Magnitudes, |Known - Calculated|", fontsize = 18)
    plt.legend()
    plt.title("Absolute Error in Star Magnitude")
    plt.savefig(f'Studentship/Graphs/FluxAccuracyHistogram_{method}.png')
    plt.show()

createFluxHistogram('mean')



