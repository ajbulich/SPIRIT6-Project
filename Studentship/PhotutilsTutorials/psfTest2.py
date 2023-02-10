import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from photutils.detection import DAOStarFinder
from astropy.visualization import LogStretch
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
# 1134,951    854, 1268    998, 988


#antiquated

sigma_psf = 2.0
image_file = 'Studentship/SPIRIT6-Observations/GoodDarks/JL-163/JL163_B_60sec_1.fts'
image_data = fits.getdata(image_file, ext=0)
image_data = image_data
w = WCS(image_file)
bkgrms = MADStdBackgroundRMS()
std = bkgrms(image_data)
daogroup = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
fitter = LevMarLSQFitter()
psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
psf_model.x_0.fixed = True
psf_model.y_0.fixed = True
sources = Table()
sources['x_mean'] = [123.5+200]
sources['y_mean'] = [58.5+200]


pos = Table(names = ['x_0', 'y_0'], data = [sources['x_mean'], sources['y_mean']])
photometry = BasicPSFPhotometry(group_maker=daogroup,
                                bkg_estimator=mmm_bkg,
                                psf_model=psf_model,
                                fitter=LevMarLSQFitter(),
                                fitshape=(9, 9))
result_tab = photometry(image=image_data, init_guesses=pos)
residual_image = photometry.get_residual_image()
print(result_tab['x_0', 'y_0', 'flux_fit'])


plt.subplot(1, 2, 1)
norm = ImageNormalize(stretch = LogStretch())
plt.imshow(image_data, aspect=1, cmap = 'gray', norm = norm, interpolation='nearest',
               origin='lower')
plt.title('Simulated data')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.subplot(1, 2, 2)
plt.imshow(residual_image, cmap='gray', aspect=1,
           interpolation='nearest',norm = norm, origin='lower')
plt.title('3.5 Sigma Threshold Subtracted')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.show()