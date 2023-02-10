import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style

#antiquated

plt.style.use(astropy_mpl_style)

from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

image_file = 'test2.fts'
fits.info(image_file)

image_data = fits.getdata(image_file, ext=0)
image_backsub = image_data - np.median(image_data)


plt.figure(1)
plt.imshow(image_data, cmap='gist_gray')
plt.colorbar()
#plt.show()

plt.figure(2)
plt.imshow(image_backsub, cmap = 'gist_gray')
plt.colorbar()
#plt.show()


from photutils.datasets import load_star_image
hdu = load_star_image()
image3 = hdu.data[500:700, 500:700].astype(float)
image2 = image3 - np.median(image3)

from photutils.detection import DAOStarFinder
from astropy.stats import mad_std

bkg_sigma = mad_std(image2)
daofind = DAOStarFinder(fwhm=4., threshold=3. * bkg_sigma)
sources = daofind(image2)
for col in sources.colnames:
    sources[col].info.format = '%.8g'
print(sources)

from photutils.aperture import aperture_photometry, CircularAperture
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=4.)
phot_table = aperture_photometry(image2, apertures)
for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'
print(phot_table)


plt.imshow(image2, cmap = 'gray_r', origin = 'lower')
apertures.plot(color = 'blue', lw = 1.5, alpha = 0.5)
plt.show()

from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
norm = ImageNormalize(stretch = SqrtStretch())
plt.imshow(image3, norm=norm, origin='lower', cmap = 'Greys_r', interpolation='nearest')
plt.show()

#using the median absolute deviation to estimate the background noise level:
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
print(mad_std(image3))
mean, median, std = sigma_clipped_stats(image3, sigma=3.0)
print((mean, median, std))