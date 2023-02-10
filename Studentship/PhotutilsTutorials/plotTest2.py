import matplotlib.pyplot as plt
from matplotlib.pyplot import*
import numpy as np
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import fits
from photutils.background import MMMBackground, MADStdBackgroundRMS

#antiquated

list60 = []
list90 = []

file = 'Studentship/PhotutilsTutorials/JL163Exposures.txt'
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
    if ExposureTime == '60sec':
        list60.append(image_file)
    elif ExposureTime == '90sec':
        list90.append(image_file)

print(list60)
print(list90)

countList60 = []
binList60 = []
order60 = []
countList90 = []
binList90 = []
order90 = []
bkgList60 = []
bkgList90 = []
rmsList60 = []
rmsList90 = []
for item in list60:
    image_data = fits.getdata(item, ext=0)
    counts3, bins3 = np.histogram(image_data, bins=16000)
    countList60.append(counts3)
    binList60.append(bins3)
    mmm_bkg = MMMBackground()
    bkg = mmm_bkg(image_data)
    bkgList60.append(bkg)
    bkgrms = MADStdBackgroundRMS()
    std = bkgrms(image_data)
    rmsList60.append(std)
    file_info = item.strip('\n').split('_')
    FilterType = file_info[1]
    order60.append(FilterType)

for item2 in list90:
    image_data = fits.getdata(item2, ext=0)
    counts2, bins2 = np.histogram(image_data, bins=16000)
    countList90.append(counts2)
    binList90.append(bins2)
    mmm_bkg = MMMBackground()
    bkg = mmm_bkg(image_data)
    bkgList90.append(bkg)
    bkgrms = MADStdBackgroundRMS()
    std = bkgrms(image_data)
    rmsList90.append(std)
    file_info = item.strip('\n').split('_')
    FilterType = file_info[1]
    order90.append(FilterType)

plt.figure(figsize = (12,8))
plt.subplot(2,2,1)
bins4 = binList60[0]
bins5 = binList60[1]
counts4 = countList60[0]
counts5 = countList60[1]
plt.hist(bins4[:-1], bins4[0:150], weights = counts4, alpha = 0.5)
plt.hist(bins5[:-1], bins5[0:150], weights = counts5, alpha = 0.5)
plt.axvline(bkgList60[0], color = 'blue')
plt.axvline(bkgList60[0] - rmsList60[0], color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkgList60[0] + rmsList60[0], color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkgList60[1], color = 'red')
plt.axvline(bkgList60[1] - rmsList60[1], color = 'red', ls = '--', alpha = 0.5)
plt.axvline(bkgList60[1] + rmsList60[1], color = 'red', ls = '--', alpha = 0.5)
plt.yscale('log')
plt.title("60 sec B Filter images")

plt.subplot(2,2,2)
bins6 = binList60[2]
bins7 = binList60[3]
counts6 = countList60[2]
counts7 = countList60[3]
plt.hist(bins6[:-1], bins6[0:150], weights = counts6, alpha = 0.5)
plt.hist(bins7[:-1], bins7[0:150], weights = counts7, alpha = 0.5)
plt.axvline(bkgList60[2], color = 'blue')
plt.axvline(bkgList60[2] - rmsList60[2], color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkgList60[2] + rmsList60[2], color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkgList60[3], color = 'red')
plt.axvline(bkgList60[3] - rmsList60[3], color = 'red', ls = '--', alpha = 0.5)
plt.axvline(bkgList60[3] + rmsList60[3], color = 'red', ls = '--', alpha = 0.5)
plt.title("60 sec V filter images")
plt.yscale('log')

plt.subplot(2,2,3)
bins8 = binList90[0]
bins9 = binList90[1]
counts8 = countList90[0]
counts9 = countList90[1]
plt.hist(bins8[:-1], bins8[0:150], weights = counts8, alpha = 0.5)
plt.hist(bins9[:-1], bins9[0:150], weights = counts9, alpha = 0.5)
plt.axvline(bkgList90[0], color = 'blue')
plt.axvline(bkgList90[0] - rmsList90[0], color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkgList90[0] + rmsList90[0], color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkgList90[1], color = 'red')
plt.axvline(bkgList90[1] - rmsList90[1], color = 'red', ls = '--', alpha = 0.5)
plt.axvline(bkgList90[1] + rmsList90[1], color = 'red', ls = '--', alpha = 0.5)
plt.title("90 sec B filter images")
plt.yscale('log')

plt.subplot(2,2,4)
bins10 = binList90[2]
bins11 = binList90[3]
counts10 = countList90[2]
counts11 = countList90[3]
plt.hist(bins10[:-1], bins10[0:150], weights = counts10, alpha = 0.5)
plt.hist(bins11[:-1], bins11[0:150], weights = counts11, alpha = 0.5)
plt.axvline(bkgList90[2], color = 'blue')
plt.axvline(bkgList90[2] - rmsList90[2], color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkgList90[2] + rmsList90[2], color = 'blue', ls = '--', alpha = 0.5)
plt.axvline(bkgList90[3], color = 'red')
plt.axvline(bkgList90[3] - rmsList90[3], color = 'red', ls = '--', alpha = 0.5)
plt.axvline(bkgList90[3] + rmsList90[3], color = 'red', ls = '--', alpha = 0.5)
plt.title("90 sec V filter images")
plt.yscale('log')


plt.savefig(f'{ObjectName}Histogram.png')
plt.show()