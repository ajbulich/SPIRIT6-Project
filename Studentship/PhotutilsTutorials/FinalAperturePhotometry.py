import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from photutils.detection import DAOStarFinder
from astropy.visualization import SqrtStretch, LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats
from astropy.wcs import WCS
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.utils import calc_total_error


def annulusPhotometry(objectName, FilterType, method, date, stackNum):

    """
    objectName - single string, e.g. JL163 or WD0830 \n
    FilterType - B or V \n
    method     - 'mean' or 'median' \n
    date       - either a single date in yyyymmdd format, or a string of multiple dates. e.g '20230125,20230127,20230128...' \n
    stackNum   - most often is 8 times the amount of dates. Just the maximum size of the stack \n
    """
    
    file2 = f'Studentship/PhotutilsTutorials/{objectName}Stars.txt'
    if method == 'median':
        file3 = f'Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_{FilterType}_{stackNum}_PSF.txt'
        file = f'Studentship/StackedImages/{date}_{objectName}Stacked_{stackNum}_{FilterType}.fits'
    else:
        file3 = f'Studentship/PhotutilsTutorials/{date}_{objectName}Stacked_{FilterType}_{stackNum}_mean_PSF.txt'
        file = f'Studentship/StackedImages/{date}_{objectName}Stacked_mean_{FilterType}_{stackNum}.fits'
    
    image_file = file
    image_data = fits.getdata(image_file, ext=0)
    w = WCS(image_file)
    sigmaList = []
    with open(file3, 'r') as f3:
        lines3 = f3.readlines()
    for line in lines3:
        var = line.strip('\n').split(' ')
        sigmaList.append(float(var[3]))
           
    sigma = np.median(sigmaList)     
    fwhm = gaussian_sigma_to_fwhm * sigma  
    i = 0
    j = 0
    with open(file2, "r") as f2:
        lines2 = f2.readlines()
    for line3 in lines2:
        j += 1
    positions2 = np.empty((j,2))
    for line2 in lines2:
        var = line2.split(' ')
        positions2[i][0] = var[0]
        positions2[i][1] = var[1]
        i+=1 
    f2.close()
    rms = MADStdBackgroundRMS()
    bkg_rms = rms(image_data)
    rms_array = bkg_rms * np.ones(np.shape(image_data), dtype = float)
    bkg = MMMBackground()
    image_bkg = bkg(image_data) * np.ones(np.shape(image_data), dtype = float)
    subtracted_image = image_data - image_bkg
    error_array = calc_total_error(subtracted_image, rms_array, 90)
    positions = w.wcs_world2pix(positions2, 1)
    apertures = CircularAperture(positions, r = 2*fwhm)
    annulus_aperture = CircularAnnulus(positions, r_in=5*fwhm, r_out=9*fwhm)
    aperstats = ApertureStats(image_data, annulus_aperture)
    bkg_mean = aperstats.mean
    phot_table = aperture_photometry(image_data, apertures, error = error_array)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    aperture_area = apertures.area_overlap(image_data)
    total_bkg = bkg_mean * aperture_area
    phot_bkgsub = phot_table['aperture_sum'] - total_bkg
    phot_error = phot_table['aperture_sum_err']
    norm = ImageNormalize(stretch = LogStretch())
    #plt.imshow(image_data, cmap='Greys', origin='lower', norm=norm,
            #interpolation='nearest')
    #apertures.plot(color='blue', lw=1.5, alpha=0.5) 
    #plt.show()
    phot_table['total_bkg'] = total_bkg
    phot_table['aperture_sum_bkgsub'] = phot_bkgsub
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    #print(phot_table)
    #write positions2[i][0] positions2[i][1] for i<j in positions2, then Exposure Filter Count. Count in phot_bkgsub[i]
    ii = 0
    if method == 'median':
        f3 = open(f"Studentship/PhotutilsTutorials/{date}_{objectName}_MedianStackedApertureResults-{FilterType}_{stackNum}.txt", "a")
        for ii in range(0,j):
            f3.write(f"{positions2[ii][0]} {positions2[ii][1]} {FilterType} {phot_bkgsub[ii]} {phot_error[ii]}\n")
        f3.close()
    else:
        f3 = open(f"Studentship/PhotutilsTutorials/{date}_{objectName}_MeanStackedApertureResults-{FilterType}_{stackNum}.txt", "a")
        for ii in range(0,j):
            f3.write(f"{positions2[ii][0]} {positions2[ii][1]} {FilterType} {phot_bkgsub[ii]} {phot_error[ii]}\n")
        f3.close()

        