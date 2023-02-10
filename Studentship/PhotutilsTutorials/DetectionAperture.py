import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import astropy_mpl_style
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
from photutils.detection import DAOStarFinder
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats
from astropy.wcs import WCS
from photutils.background import MMMBackground, MADStdBackgroundRMS

#THIS IS ANTIQUATED

image_file = 'Studentship/SPIRIT6-Observations/GoodDarks/JL-163/JL163_B_60sec_1.fts'
fits.info(image_file)
file_info2 = image_file.strip('\n').split('/')
file_info = file_info2[4]
ObjectName = file_info[0]
FilterType = file_info[1]
ExposureTime = file_info[2]

image_data = fits.getdata(image_file, ext=0)
w = WCS(image_file)
bkgrms = MADStdBackgroundRMS()
std = bkgrms(image_data)
mmm_bkg = MMMBackground()
daofind = DAOStarFinder(threshold = ( mmm_bkg(image_data)+2.5*std), fwhm = 2)
sources = daofind(image_data - mmm_bkg(image_data))
positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=7.)
norm = ImageNormalize(stretch = SqrtStretch())
annulus_aperture = CircularAnnulus(positions, r_in=10, r_out=15)
aperstats = ApertureStats(image_data, annulus_aperture)
bkg_mean = aperstats.mean
phot_table = aperture_photometry(image_data, apertures)
for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output
aperture_area = apertures.area_overlap(image_data)
total_bkg = bkg_mean * aperture_area
phot_bkgsub = phot_table['aperture_sum'] - total_bkg
phot_table['total_bkg'] = total_bkg
phot_table['aperture_sum_bkgsub'] = phot_bkgsub
for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output

ascii.write(phot_table[...], format = 'csv', output = f'{ObjectName}-{ExposureTime}Aperture.csv')
print(phot_table)