from astropy.visualization import simple_norm
import matplotlib.pyplot as plt
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.datasets import make_100gaussians_image
from astropy.io import fits
from hourstodegrees import*

#used to create example image for photometry.

file1 = 'TPhe_V_45sec_1.fts'
data = fits.getdata(file1,ext=0)
positions = [(1219.6, 1026.58), (794.168, 1588.94), (1266.95, 1272.93)]
aperture = CircularAperture(positions, r=8)
ap_patches = aperture.plot(color='white', lw=2,
                           label='Photometry aperture')
annulus_aperture = CircularAnnulus(positions, r_in=15, r_out=23)

norm = simple_norm(data, 'sqrt', percent=99)
plt.imshow(data, norm=norm, interpolation='nearest')
plt.xlim(1100, 1400)
plt.ylim(1000, 1300)


ann_patches = annulus_aperture.plot(color='red', lw=2,
                                    label='Background annulus')
handles = (ap_patches[0], ann_patches[0])
plt.legend(loc=(0.3, 0.2), facecolor='#458989', labelcolor='white',
           handles=handles, prop={'weight': 'bold', 'size': 11})
plt.show()