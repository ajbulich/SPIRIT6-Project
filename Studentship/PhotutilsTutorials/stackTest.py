
from distutils.log import Log
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from StackWeight import StackWeighting
from photutils.background import MMMBackground, MADStdBackgroundRMS

#antiquated

file = 'Studentship/PhotutilsTutorials/JL163ReprojectedImages-V.txt'

RMSArr = StackWeighting('JL163', 'V')
with open(file, 'r') as f:
    lines = f.readlines()

data_arr = []
for line in lines:
    image_file = line.strip('\n')
    image_data = fits.getdata(image_file, ext=0)
    data_arr.append(image_data)

minRMS = float(np.amin(RMSArr))
newArr = (1 / minRMS) * np.array( RMSArr)
avg = float(np.average(newArr))
newArr2 = []

for item in newArr:
    newArr2.append(avg / item)

print(newArr2)

newData = []

for i in range(0,len(newArr2)):
    newData.append(newArr2[i] * data_arr[i])

d = np.stack(newData)
final = np.median(d, axis=0)

norm = ImageNormalize(stretch = LogStretch())
plt.imshow(final, cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
plt.show()

plt.imshow(data_arr[1], cmap='Greys', origin='lower', norm=norm, interpolation='nearest')
plt.show()

bkgrms = MADStdBackgroundRMS()
print(bkgrms(final))

RMSList = []
for item in data_arr:
    RMSList.append(bkgrms(item))

print(RMSList)

file = 'Studentship/StackedImages/JL163Stacked_8_V.fits'
image_data = fits.getdata(file, ext=0)
print(bkgrms(image_data))


#compute Background rms for each image beforehand, and then do it afterward on the stacked image too. Compare and hope its lower.
#generate a histogram of unstacked vs stacked.