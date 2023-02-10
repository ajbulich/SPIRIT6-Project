import numpy as np
import matplotlib.pyplot as plt

#ANTIQUATED


file = 'Studentship/PhotutilsTutorials/JL163ApertureResults.txt'
B60FluxList1 = []
B90FluxList1 = []
V60FluxList1 = []
V90FluxList1 = []
with open(file, "r") as f:
    lines = f.readlines()
for line in lines:
    var = line.strip('\n').split(' ')
    if((var[2] == '60sec') and (var[3] == 'B')):
        B60FluxList1.append(float(var[4]))
    elif((var[2] == '90sec') and (var[3] == 'B')):
        B90FluxList1.append(float(var[4]))
    elif((var[2] == '60sec') and (var[3] == 'V')):
        V60FluxList1.append(float(var[4]))
    elif((var[2] == '90sec') and (var[3] == 'V')):
        V90FluxList1.append(float(var[4]))
f.close()

file2 = 'Studentship/PhotutilsTutorials/JL163PSF.txt'
B60FluxList2 = []
B60SigmaList2 = []
B90FluxList2 = []
B90SigmaList2 = []
V60FluxList2 = []
V60SigmaList2 = []
V90FluxList2 = []
V90SigmaList2 = []
with open(file2, "r") as f2:
    lines2 = f2.readlines()
for line2 in lines2:
    var2 = (line2.strip('\n').split(' '))
    if((var2[2] == '60sec') and (var2[3] == 'B')):
        B60FluxList2.append(float(var2[4]))
        B60SigmaList2.append(float(var2[5]))
    elif((var2[2] == '90sec') and (var2[3] == 'B')):
        B90FluxList2.append(float(var2[4]))
        B90SigmaList2.append(float(var2[5]))
    elif((var2[2] == '60sec') and (var2[3] == 'V')):
        V60FluxList2.append(float(var2[4]))
        V60SigmaList2.append(float(var2[5]))
    elif((var2[2] == '90sec') and (var2[3] == 'V')):
        V90FluxList2.append(float(var2[4]))
        V90SigmaList2.append(float(var2[5]))
f2.close()
x = np.linspace(0,65000,10)
y = x

plt.subplot(2,2,1)
plt.axhline(0)
plt.scatter(B60FluxList1, np.log(B60FluxList2)-np.log(B60FluxList1))
plt.xlabel('Aperture Flux')
plt.ylabel('log(PSF) - log(Aperture)')
plt.title("60 sec B Filter")

x = np.linspace(0,100000,10)
y = x
plt.subplot(2,2,2)
plt.axhline(0)
plt.scatter(B90FluxList1, np.log(B90FluxList2)-np.log(B90FluxList1))
plt.xlabel('Aperture Flux')
plt.ylabel('log(PSF) - log(Aperture)')
plt.title("90 sec B Filter")

x = np.linspace(0,100000,10)
y = x
plt.subplot(2,2,3)
plt.axhline(0)
plt.scatter(V60FluxList1, np.log(V60FluxList2)-np.log(V60FluxList1))
plt.xlabel('Aperture Flux')
plt.ylabel('log(PSF) - log(Aperture)')
plt.title("60 sec V Filter")

x = np.linspace(0,150000,10)
y = x
plt.subplot(2,2,4)
plt.axhline(0)
plt.scatter(V90FluxList1, np.log(V90FluxList2) - np.log(V90FluxList1))
plt.xlabel('Aperture Flux')
plt.ylabel('log(PSF) - log(Aperture)')
plt.title("90 sec V Filter")

plt.show()


plt.subplot(2,2,1)
counts2, bins2 = np.histogram(B60SigmaList2, bins = 20)
plt.hist(bins2[:-1], bins2, weights = counts2, alpha = 0.5)
plt.title("B Filter 60 second")
plt.xlabel("Sigma (pixels)")
plt.ylabel('counts')

plt.subplot(2,2,2)
counts2, bins2 = np.histogram(B60SigmaList2, bins = 20)
plt.hist(bins2[:-1], bins2, weights = counts2, alpha = 0.5)
plt.title("B Filter 90 second")
plt.xlabel("Sigma (pixels)")
plt.ylabel('counts')

plt.subplot(2,2,3)
counts2, bins2 = np.histogram(B60SigmaList2, bins = 20)
plt.hist(bins2[:-1], bins2, weights = counts2, alpha = 0.5)
plt.title("V Filter 60 second")
plt.xlabel("Sigma (pixels)")
plt.ylabel('counts')

plt.subplot(2,2,4)
counts2, bins2 = np.histogram(B60SigmaList2, bins = 20)
plt.hist(bins2[:-1], bins2, weights = counts2, alpha = 0.5)
plt.title("V Filter 90 second")
plt.xlabel("Sigma (pixels)")
plt.ylabel('counts')

plt.show()





