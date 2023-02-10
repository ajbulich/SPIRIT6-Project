
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


#ANTIQUATED

def AnnulusPhotometry():
    file = 'Studentship/PhotutilsTutorials/JL163StackedImages.txt'
    file2 = 'Studentship/PhotutilsTutorials/JL163Stars.txt'
    file3 = 'Studentship/PhotutilsTutorials/JL163StackedPSF'
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
        #mean, median, std = sigma_clipped_stats(image_data, sigma = 3.0)
        #daofind = DAOStarFinder(fwhm = 3.0, threshold = 5.*std)
        #sources = daofind(image_data - median)
        #positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
           
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
    
        #positions = w.wcs_world2pix(97.285125, -71.494161, 1)   #comment this out if you want all the stars in the image
        positions = w.wcs_world2pix(positions2, 1)
        apertures = CircularAperture(positions, r=7.)
        norm = ImageNormalize(stretch = SqrtStretch())
        plt.imshow(image_data, cmap='Greys', origin='lower', norm=norm,
            interpolation='nearest')
        apertures.plot(color='blue', lw=1.5, alpha=0.5) 
        #plt.show()
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
        #print(phot_table)
        #write positions2[i][0] positions2[i][1] for i<j in positions2, then Exposure Filter Count. Count in phot_bkgsub[i]
        ii = 0
        f3 = open(f"{ObjectName}ApertureResults.txt", "a")
        for ii in range(0,j):
            f3.write(f"{positions2[ii][0]} {positions2[ii][1]} {ExposureTime} {FilterType} {phot_bkgsub[ii]}\n")
        f3.close()
        

        
    f.close()

    
    

AnnulusPhotometry()