from astropy.io import fits
from reproject import reproject_interp, reproject_adaptive, reproject_exact
from astropy.wcs import WCS
import matplotlib.pyplot as plt

#antiquated

hdu1 = fits.open('Studentship/PhotutilsTutorials/WD-0830-535_B_90sec_1.fts')[0]
hdu2 = fits.open('Studentship/PhotutilsTutorials/WD-0830-535_B_90sec_2.fts')[0]
array, footprint = reproject_exact(hdu2, hdu1.header)


fits.writeto('Studentship/PhotutilsTutorials/WD-0830-535_B_90sec_stacked3.fits', array, hdu1.header, overwrite=True)
