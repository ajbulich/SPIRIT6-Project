import numpy as np
import matplotlib.pyplot as plt

#antiquated


file1 = 'Studentship/PhotutilsTutorials/JL163StackedPSF.txt'
file2 = 'Studentship/PhotutilsTutorials/WD0830StackedPSF.txt'

with open(file1, 'r') as f:
    lines = f.readlines()

JL163Stars = 7
i = 0
originalFlux = []
originalSigma = []
Flux2 = []
Sigma2 = []
Flux3 = []
Sigma3 = []
Flux4 = []
Sigma4 = []
Flux5 = []
Sigma5 = []
Flux6 = []
Sigma6 = []
Flux7 = []
Sigma7 = []
Flux7 = []
Sigma7 = []
Flux8 = []
Sigma8 = []
for line in lines:
    var = line.strip('\n').split(' ')
    if (np.floor(i/7) == 0):
        originalFlux.append(float(var[2]))
        originalSigma.append(float(var[3]))
    elif (np.floor(i/7) == 1):
        Flux2.append(float(var[2]))
        Sigma2.append(float(var[3]))
    elif (np.floor(i/7) == 2):
        Flux3.append(float(var[2]))
        Sigma3.append(float(var[3]))
    elif (np.floor(i/7) == 3):
        Flux4.append(float(var[2]))
        Sigma4.append(float(var[3]))
    elif (np.floor(i/7) == 4):
        Flux5.append(float(var[2]))
        Sigma5.append(float(var[3]))
    elif (np.floor(i/7) == 5):
        Flux6.append(float(var[2]))
        Sigma6.append(float(var[3]))
    elif (np.floor(i/7) == 6):
        Flux7.append(float(var[2]))
        Sigma7.append(float(var[3]))
    elif (np.floor(i/7) == 7):
        Flux8.append(float(var[2]))
        Sigma8.append(float(var[3]))
    i += 1


meanArr = np.mean([originalFlux, Flux2, Flux3, Flux4, Flux5, Flux6, Flux7, Flux8], axis = 0)
meanArr2 = np.mean([originalSigma, Sigma2, Sigma3, Sigma4, Sigma5, Sigma6, Sigma7, Sigma8],axis = 0)
plt.plot(originalFlux, 'o')
plt.plot(Flux2, 'o')
plt.plot(Flux3, 'o')
plt.plot(Flux4, 'o')
plt.plot(Flux5, 'o')
plt.plot(Flux6, 'o')
plt.plot(Flux7, 'o')
plt.plot(Flux8, 'o')
plt.plot(meanArr)
plt.legend(['Original', '2 images', '3 images', '4 images', '5 images', '6 images', '7 images', '8 images', 'mean'])
plt.ylabel('Flux (ADU)')
plt.xlabel('Stars')
plt.title("JL 163 Flux for N stacked images")

plt.show()

plt.plot(originalSigma, 'o')
plt.plot(Sigma2, 'o')
plt.plot(Sigma3, 'o')
plt.plot(Sigma4, 'o')
plt.plot(Sigma5, 'o')
plt.plot(Sigma6, 'o')
plt.plot(Sigma7, 'o')
plt.plot(Sigma8, 'o')
plt.plot(meanArr2)
plt.legend(['Original', '2 images', '3 images', '4 images', '5 images', '6 images', '7 images', '8 images', 'mean'])
plt.ylabel('Sigma (Pixels)')
plt.xlabel('Stars')
plt.title("JL163 Sigma for N stacked images")

plt.show()







with open(file2, 'r') as f2:
    lines2 = f2.readlines()

JL163Stars = 7
i = 0
originalFlux = []
originalSigma = []
Flux2 = []
Sigma2 = []
Flux3 = []
Sigma3 = []
Flux4 = []
Sigma4 = []
Flux5 = []
Sigma5 = []
Flux6 = []
Sigma6 = []
Flux7 = []
Sigma7 = []
Flux7 = []
Sigma7 = []
Flux8 = []
Sigma8 = []
for line2 in lines2:
    var = line2.strip('\n').split(' ')
    if (np.floor(i/12) == 0):
        originalFlux.append(float(var[2]))
        originalSigma.append(float(var[3]))
    elif (np.floor(i/12) == 1):
        Flux2.append(float(var[2]))
        Sigma2.append(float(var[3]))
    elif (np.floor(i/12) == 2):
        Flux3.append(float(var[2]))
        Sigma3.append(float(var[3]))
    elif (np.floor(i/12) == 3):
        Flux4.append(float(var[2]))
        Sigma4.append(float(var[3]))
    elif (np.floor(i/12) == 4):
        Flux5.append(float(var[2]))
        Sigma5.append(float(var[3]))
    elif (np.floor(i/12) == 5):
        Flux6.append(float(var[2]))
        Sigma6.append(float(var[3]))
    elif (np.floor(i/12) == 6):
        Flux7.append(float(var[2]))
        Sigma7.append(float(var[3]))
    elif (np.floor(i/12) == 7):
        Flux8.append(float(var[2]))
        Sigma8.append(float(var[3]))
    i += 1




meanArr = np.mean([originalFlux, Flux2, Flux3, Flux4, Flux5, Flux6, Flux7, Flux8], axis = 0)
meanArr2 = np.mean([originalSigma, Sigma2, Sigma3, Sigma4, Sigma5, Sigma6, Sigma7, Sigma8],axis = 0)
plt.plot(originalFlux, 'o')
plt.plot(Flux2, 'o')
plt.plot(Flux3, 'o')
plt.plot(Flux4, 'o')
plt.plot(Flux5, 'o')
plt.plot(Flux6, 'o')
plt.plot(Flux7, 'o')
plt.plot(Flux8, 'o')
plt.plot(meanArr)
plt.legend(['Original', '2 images', '3 images', '4 images', '5 images', '6 images', '7 images', '8 images', 'mean'])
plt.ylabel('Flux (ADU)')
plt.xlabel('Stars')
plt.title("WD0830 Flux for N Stacked Images")

plt.show()

plt.plot(originalSigma, 'o')
plt.plot(Sigma2, 'o')
plt.plot(Sigma3, 'o')
plt.plot(Sigma4, 'o')
plt.plot(Sigma5, 'o')
plt.plot(Sigma6, 'o')
plt.plot(Sigma7, 'o')
plt.plot(Sigma8, 'o')
plt.plot(meanArr2)
plt.legend(['Original', '2 images', '3 images', '4 images', '5 images', '6 images', '7 images', '8 images', 'mean'])
plt.ylabel('Sigma (Pixels)')
plt.xlabel('Stars')
plt.title("WD0830 Sigma for N Stacked Images")

plt.show()