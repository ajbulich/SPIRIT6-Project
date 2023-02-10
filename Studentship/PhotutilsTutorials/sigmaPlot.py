import numpy as np
import matplotlib.pyplot as plt

#antiquated

file = 'Studentship/PhotutilsTutorials/TPhe_45sec_stacked_DAOPSF.txt'
sigmaList = []
fluxList = []
with open(file, "r") as f:
    lines = f.readlines()
for line in lines:
    var = line.strip('\n').split(' ')
    if (float(var[2]) <= 2000000) and (float(var[2]) >= 0):
        fluxList.append(float(var[2]))
        sigmaList.append(float(var[3]))

print(fluxList)

plt.figure(figsize = (12,8))
plt.subplot(1,3,1)
counts2, bins2 = np.histogram(sigmaList, bins = 40)
plt.hist(bins2[:-1], bins2, weights = counts2, alpha = 0.5)
plt.title("Sigma Histogram")
plt.xlabel("Sigma (pixels)")
plt.ylabel('counts')

plt.subplot(1,3,2)
plt.scatter(sigmaList, np.log(fluxList))
plt.title("log(PSF Flux) vs Sigma")
plt.xlabel('Sigma (pixels)')
plt.ylabel('log(PSF Flux)')

plt.subplot(1,3,3)
counts2, bins2 = np.histogram(fluxList, bins = 800)
plt.hist(bins2[:-1], bins2[0:40], weights = counts2, alpha = 0.5)
plt.title("Flux Histogram")
plt.xlabel("Flux")
plt.ylabel('counts')

plt.show()