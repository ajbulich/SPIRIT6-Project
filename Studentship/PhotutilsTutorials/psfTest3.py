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

#antiquated

def PSFPhotometry():
    file = 'Studentship/PhotutilsTutorials/JL163StackedImages-V.txt'
    file2 = 'Studentship/PhotutilsTutorials/JL163Stars.txt'
    plt.style.use(astropy_mpl_style)
    with open(file, "r") as f:
        lines = f.readlines()
    for line in lines:
        image_file = line.strip('\n')
        fits.info(image_file)
        file_info = line.strip('\n').split('_')
        ObjectName = file_info[0]
        FilterType = file_info[1]
        ExposureTime = file_info[2]
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
        gaussian_PRF.sigma.value = 1.1
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
        f3 = open(f"{ObjectName}PSFTesting.txt", "a")
        for ii in range(0,j):
            f3.write(f"{positions2[ii][0]} {positions2[ii][1]} {ExposureTime} {FilterType} {result_tab['flux_fit'][ii]} {result_tab['sigma_fit'][ii]}\n")
        f3.close()
        

    f.close()

PSFPhotometry()
    