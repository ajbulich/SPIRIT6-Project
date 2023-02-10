import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from photutils.detection import DAOStarFinder
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats
from astropy.wcs import WCS
from astropy.table import Table
from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.psf import IterativelySubtractedPSFPhotometry
from astropy.stats import gaussian_sigma_to_fwhm
from photutils.psf import BasicPSFPhotometry

def SingleImagePSF(image_file, objectName):
    """
    Designed to be used in conjuction with loop through the 90secExposures file for a given object/date.\n
    \n
    objectName - single object name string, e.g 'JL163'
    image_file - exact path to specific image
    Returns - \n
        - list() type that contains all sigma values for that object/date combo\n
    """

    file2 = f'Studentship/PhotutilsTutorials/{objectName}Stars.txt'
    image_data = fits.getdata(image_file, ext=0)
    w = WCS(image_file)
    j = 0
    i = 0
    with open(file2, 'r') as f:
        lines = f.readlines()
    for line in lines:
        j += 1
    positions2 = np.empty((j,2))
    for line2 in lines:
        var = line2.split(' ')
        positions2[i][0] = var[0]
        positions2[i][1] = var[1]
        i+=1 
    f.close()

    positions = w.wcs_world2pix(positions2, 1)
    print(positions)
    bkgrms = MADStdBackgroundRMS()
    std = bkgrms(image_data)
    daogroup = DAOGroup(crit_separation = 8.0)
    mmm_bkg = MMMBackground()
    fitter = LevMarLSQFitter()
    gaussian_PRF = IntegratedGaussianPRF()
    gaussian_PRF.sigma.fixed = False
    gaussian_PRF.sigma.value = 2.2
    gaussian_PRF.x_0.fixed = True
    gaussian_PRF.y_0.fixed = True
    xlist = []
    ylist = []
    for item in positions:
        xlist.append(item[0])
        ylist.append(item[1])
    sources = Table()
    sources['x_mean'] = xlist
    sources['y_mean'] = ylist
    pos = Table(names = ['x_0', 'y_0'], data = [sources['x_mean'], sources['y_mean']])
    photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=mmm_bkg,
                                    psf_model=gaussian_PRF,
                                    fitter=LevMarLSQFitter(),
                                    fitshape=(9, 9))
    result_tab = photometry(image=image_data, init_guesses=pos)
    sigma_arr = result_tab['sigma_fit']
    return sigma_arr



def PSFStackedPhotometry(numImage, FilterType, ObjectName, method, date, stackNum):
    if numImage == 1:
        if method == 'mean':
            image_file = f'Studentship/StackedImages/{date}_{ObjectName}Stacked_mean_{FilterType}_{stackNum}.fits'
        else:
            image_file = f'Studentship/StackedImages/{date}_{ObjectName}Stacked_{stackNum}_{FilterType}.fits'
        file2 = f'Studentship/PhotutilsTutorials/{ObjectName}Stars.txt'
        image_data = fits.getdata(image_file, ext = 0)
        w = WCS(image_file)
        i = 0
        j = 0
        with open(file2, 'r') as f:
            lines = f.readlines()
        for line in lines:
            j += 1
        positions2 = np.empty((j,2))
        for line2 in lines:
            var = line2.split(' ')
            positions2[i][0] = var[0]
            positions2[i][1] = var[1]
            i+=1 
        f.close()

        positions = w.wcs_world2pix(positions2, 1)
        bkgrms = MADStdBackgroundRMS()
        std = bkgrms(image_data)
        daogroup = DAOGroup(crit_separation = 8.0)
        mmm_bkg = MMMBackground()
        fitter = LevMarLSQFitter()
        gaussian_PRF = IntegratedGaussianPRF()
        gaussian_PRF.sigma.fixed = False
        gaussian_PRF.sigma.value = 2.2
        gaussian_PRF.x_0.fixed = True
        gaussian_PRF.y_0.fixed = True
        xlist = []
        ylist = []
        for item in positions:
            xlist.append(item[0])
            ylist.append(item[1])
        sources = Table()
        sources['x_mean'] = xlist
        sources['y_mean'] = ylist
        pos = Table(names = ['x_0', 'y_0'], data = [sources['x_mean'], sources['y_mean']])
        print(pos)
        photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=mmm_bkg,
                                    psf_model=gaussian_PRF,
                                    fitter=LevMarLSQFitter(),
                                    fitshape=(9, 9))
        result_tab = photometry(image=image_data, init_guesses=pos)
        print(result_tab['flux_fit', 'sigma_fit'])
        ii = 0
        if method == 'median':
            f3 = open(f"Studentship/PhotutilsTutorials/{date}_{ObjectName}Stacked_{FilterType}_{stackNum}_PSF.txt", "a")
            for ii in range(0,j):
                f3.write(f"{positions2[ii][0]} {positions2[ii][1]} {result_tab['flux_fit'][ii]} {result_tab['sigma_fit'][ii]}\n")
            f3.close()
        else:
            f3 = open(f"Studentship/PhotutilsTutorials/{date}_{ObjectName}Stacked_{FilterType}_{stackNum}_mean_PSF.txt", "a")
            for ii in range(0,j):
                f3.write(f"{positions2[ii][0]} {positions2[ii][1]} {result_tab['flux_fit'][ii]} {result_tab['sigma_fit'][ii]}\n")
            f3.close()


    elif numImage > 1:
        file = f'Studentship/StackedImages/{date}_{ObjectName}StackedImages-{FilterType}.txt'
        file2 = f'Studentship/PhotutilsTutorials/{ObjectName}Stars.txt'
        plt.style.use(astropy_mpl_style)
        with open(file, "r") as f:
            lines = f.readlines()
        for line in lines:
            image_file = line.strip('\n')
            fits.info(image_file)
            image_data = fits.getdata(image_file, ext=0)
            w = WCS(image_file)
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
            positions = w.wcs_world2pix(positions2, 1)
            bkgrms = MADStdBackgroundRMS()
            std = bkgrms(image_data)
            daogroup = DAOGroup(crit_separation = 8.0)
            mmm_bkg = MMMBackground()
            fitter = LevMarLSQFitter()
            gaussian_PRF = IntegratedGaussianPRF()
            gaussian_PRF.sigma.fixed = False
            gaussian_PRF.sigma.value = 2.2
            gaussian_PRF.x_0.fixed = True
            gaussian_PRF.y_0.fixed = True
            xlist = []
            ylist = []
            for item in positions:
                xlist.append(item[0])
                ylist.append(item[1])
            sources = Table()
            sources['x_mean'] = xlist
            sources['y_mean'] = ylist
            pos = Table(names = ['x_0', 'y_0'], data = [sources['x_mean'], sources['y_mean']])
            photometry = BasicPSFPhotometry(group_maker=daogroup,
                                    bkg_estimator=mmm_bkg,
                                    psf_model=gaussian_PRF,
                                    fitter=LevMarLSQFitter(),
                                    fitshape=(9, 9))
            result_tab = photometry(image=image_data, init_guesses=pos)
            print(result_tab['flux_fit', 'sigma_fit'])
            ii = 0
            f3 = open(f"Studentship/PhotutilsTutorials/{date}_{ObjectName}StackedPSF_{FilterType}.txt", "a")
            for ii in range(0,j):
                f3.write(f"{positions2[ii][0]} {positions2[ii][1]} {result_tab['flux_fit'][ii]} {result_tab['sigma_fit'][ii]}\n")
            f3.close()
            

        f.close()

    