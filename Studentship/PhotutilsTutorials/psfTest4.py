
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

gaussian_PRF = IntegratedGaussianPRF()
gaussian_PRF.sigma.fixed = False
gaussian_PRF.sigma.value = 1.1

image_file = 'Studentship/PhotutilsTutorials/JL163Stacked_8_B.fits'
image_data2 = fits.getdata(image_file, ext=0)
image_data = image_data2

w = WCS(image_file)
fits.info(image_file)
bkgrms = MADStdBackgroundRMS()
std = bkgrms(image_data)
mmm_bkg = MMMBackground()
iraffind = DAOStarFinder(threshold = ( mmm_bkg(image_data)+2.5*std), fwhm = 2)
daogroup = DAOGroup(crit_separation = 8.0)
fitter = LevMarLSQFitter()
sources = iraffind(image_data-mmm_bkg(image_data))
photometry = IterativelySubtractedPSFPhotometry(finder=iraffind,group_maker=daogroup,bkg_estimator=mmm_bkg ,psf_model=gaussian_PRF,fitter=LevMarLSQFitter(),niters=2, fitshape=(11, 11))
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



print(result_tab['x_fit', 'y_fit', 'flux_fit', 'sigma_fit'])
rows = (np.shape(result_tab))[0]

f3 = open(f"Studentship/PhotutilsTutorials/WD-0830-535_B_90sec_stacked3_DAOPSF.txt", "a")
for ii in range(0,rows):
    f3.write(f"{result_tab[ii]['x_fit']} {result_tab[ii]['y_fit']} {result_tab[ii]['flux_fit']} {result_tab[ii]['sigma_fit']}\n")
f3.close()

