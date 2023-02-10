import matplotlib.pyplot as plt
from matplotlib.pyplot import*
import numpy as np
from astropy.visualization import astropy_mpl_style
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from photutils.detection import DAOStarFinder
from astropy.visualization import SqrtStretch
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats
from astropy.wcs import WCS
from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from photutils.psf import IterativelySubtractedPSFPhotometry
from astropy.stats import gaussian_sigma_to_fwhm

#antiquated

sigma_psf = 2.0
image_file = 'Studentship/SPIRIT6-Observations/GoodDarks/JL-163/JL163_B_60sec_1.fts'
image_data2 = fits.getdata(image_file, ext=0)
image_data = image_data2
print(np.shape(image_data))
w = WCS(image_file)
bkgrms = MADStdBackgroundRMS()
std = bkgrms(image_data)
mmm_bkg = MMMBackground()
iraffind = DAOStarFinder(threshold = ( mmm_bkg(image_data)), fwhm = sigma_psf * gaussian_sigma_to_fwhm, roundhi=5.0, roundlo=-5.0, sharphi=2.0, sharplo=0.0)
daogroup = DAOGroup(2.0 * sigma_psf * gaussian_sigma_to_fwhm)
fitter = LevMarLSQFitter()
sources = iraffind(image_data-mmm_bkg(image_data))
print(type(sources['xcentroid', 'ycentroid'][0]))
psf_model = IntegratedGaussianPRF(sigma=sigma_psf)
photometry = IterativelySubtractedPSFPhotometry(finder=iraffind,group_maker=daogroup,bkg_estimator=mmm_bkg ,psf_model=psf_model,fitter=LevMarLSQFitter(),niters=1, fitshape=(11, 11))
result_tab = photometry(image=image_data)
residual_image = photometry.get_residual_image()

plt.subplot(1, 2, 1)
norm = ImageNormalize(stretch = LogStretch())

plt.imshow(image_data, cmap='gray', aspect=1,norm = norm, interpolation='nearest', origin='lower')
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=6.)
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.title('Simulated data')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.subplot(1, 2, 2)
plt.imshow(residual_image, cmap='gray', aspect=1, interpolation='nearest', origin='lower')
plt.title('Residual Image')
plt.colorbar(orientation='horizontal', fraction=0.046, norm = norm, pad=0.04)

plt.show()


print(result_tab['x_fit', 'y_fit', 'flux_fit'])

