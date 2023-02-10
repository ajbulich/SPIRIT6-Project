import matplotlib.pyplot as plt
from matplotlib.pyplot import*
import numpy as np
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from photutils.background import MMMBackground, MADStdBackgroundRMS

#antiquated

image_file = 'Studentship/StackedImages/20230116_JL163Stacked_mean_B_8.fits'
image_data = fits.getdata(image_file, ext=0)[100:700,480:1230]

counts3, bins3 = np.histogram(image_data, bins=16000)
mmm_bkg = MMMBackground()
bkg1 = mmm_bkg(image_data)
bkgrms = MADStdBackgroundRMS()
rms1 = bkgrms(image_data)

image_file2 = 'Studentship/StackedImages/20230115_JL163Stacked_mean_B_8.fits'
image_data2 = fits.getdata(image_file2, ext=0)[100:700,480:1230]
counts4, bins4 = np.histogram(image_data2, bins=16000)
mmm_bkg2 = MMMBackground()
bkg2 = mmm_bkg2(image_data2)
bkgrms2 = MADStdBackgroundRMS()
rms2 = bkgrms2(image_data2)

image_file3 = 'Studentship/StackedImages/20221214_JL163Stacked_mean_B_8.fits'
image_data3 = fits.getdata(image_file3, ext=0)[100:700,480:1230]
counts5, bins5 = np.histogram(image_data3, bins=16000)
mmm_bkg3 = MMMBackground()
bkg3 = mmm_bkg3(image_data3)
bkgrms3 = MADStdBackgroundRMS()
rms3 = bkgrms3(image_data3)


plt.figure(figsize = (16,8))
plt.subplot(1,3,1)   #COMPARE 1 to 2
plt.hist(bins3[:-1], bins3[0:200], weights = counts3, alpha = 0.5)
plt.hist(bins4[:-1], bins4[0:200], weights = counts4, alpha = 0.5)
plt.axvline(bkg1, color = 'blue')
plt.axvline(bkg1 - rms1, color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkg1 + rms1, color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkg2, color = 'red')
plt.axvline(bkg2 - rms2, color = 'red', ls = '--', alpha = 0.5)
plt.axvline(bkg2 + rms2, color = 'red', ls = '--', alpha = 0.5)
plt.yscale('log')
plt.title("1 compared to 2")

plt.subplot(1,3,2)   #COMPARE 1 to 2
plt.hist(bins4[:-1], bins4[0:200], weights = counts4, alpha = 0.5)
plt.hist(bins5[:-1], bins5[0:200], weights = counts5, alpha = 0.5)
plt.axvline(bkg2, color = 'blue')
plt.axvline(bkg2 - rms2, color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkg2 + rms2, color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkg3, color = 'red')
plt.axvline(bkg3 - rms3, color = 'red', ls = '--', alpha = 0.5)
plt.axvline(bkg3 + rms3, color = 'red', ls = '--', alpha = 0.5)
plt.yscale('log')
plt.title("2 compared to 3")

plt.subplot(1,3,3)   #COMPARE 1 to 2
plt.hist(bins3[:-1], bins3[0:200], weights = counts3, alpha = 0.5)
plt.hist(bins5[:-1], bins5[0:200], weights = counts5, alpha = 0.5)
plt.axvline(bkg1, color = 'blue')
plt.axvline(bkg1 - rms1, color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkg1 + rms1, color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkg3, color = 'red')
plt.axvline(bkg3 - rms3, color = 'red', ls = '--', alpha = 0.5)
plt.axvline(bkg3 + rms3, color = 'red', ls = '--', alpha = 0.5)
plt.yscale('log')
plt.title("1 compared to 3")



plt.show()


plt.hist(bins4[:-1], bins4[0:200], weights = counts4, alpha = 0.5)
plt.hist(bins5[:-1], bins5[0:200], weights = counts5, alpha = 0.5)
plt.axvline(bkg2, color = 'blue')
plt.axvline(bkg2 - rms2, color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkg2 + rms2, color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkg3, color = 'red')
plt.axvline(bkg3 - rms3, color = 'red', ls = '--', alpha = 0.5)
plt.axvline(bkg3 + rms3, color = 'red', ls = '--', alpha = 0.5)
plt.yscale('log')
plt.title("2 compared to 3")


plt.show()