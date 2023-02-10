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
image_file = 'Studentship/SPIRIT6-Observations/GoodDarks/JL-82/JL-82_B_60sec_1.fts'
image_data2 = fits.getdata(image_file, ext=0)
image_data = image_data2
print(np.shape(image_data))
w = WCS(image_file)
bkgrms = MADStdBackgroundRMS()
std = bkgrms(image_data)
threshold = 3.5*std
reducedImage = image_data - threshold


plt.subplot(1, 2, 1)
norm = ImageNormalize(stretch = LogStretch())
plt.imshow(image_data, cmap='viridis', aspect=1,norm = norm, interpolation='nearest', origin='lower')
plt.title('Normal Image')
plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
plt.subplot(1, 2, 2)
plt.imshow(reducedImage, cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.title('Reduced')
plt.colorbar(orientation='horizontal', fraction=0.046, norm = norm, pad=0.04)
plt.show()



