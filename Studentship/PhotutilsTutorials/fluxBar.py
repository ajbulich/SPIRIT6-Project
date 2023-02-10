import numpy as np
import matplotlib.pyplot as plt
import math

#antiquated

file = 'Studentship/PhotutilsTutorials/WD-0830-535PSFTesting.txt'
file2 = 'Studentship/PhotutilsTutorials/WD0830Stars.txt'
with open(file2, "r") as f2:
    lines = f2.readlines()
    j = 0
    for line in lines:
        j += 1

print(j)
numStars = j
Flux1 = []
Sigma1 = []
Flux2 = []
Sigma2 = []
Flux3 = []
Sigma3 = []
Flux4 = []
Sigma4 = []
with open(file, "r") as f:
    lines = f.readlines()
    i = 0
    for line in lines:
        var = line.strip('\n').split(' ')
        flux = float(var[4])
        sigma = float(var[5])
        if math.floor(i/j) == 0:
            Flux1.append(flux)
            Sigma1.append(sigma)
        elif math.floor(i/j) == 1:
            Flux2.append(flux)
            Sigma2.append(sigma)
        elif math.floor(i/j) == 2:
            Flux3.append(flux)
            Sigma3.append(sigma)
        elif(math.floor(i/j) == 3):
            Flux4.append(flux)
            Sigma4.append(sigma)
        i += 1

Result1 = []
Result2 = []
Result3 = []
for i in range(0,numStars):
    Result1.append(Flux1[i]/Flux4[i])
    Result2.append(Flux2[i]/Flux4[i])
    Result3.append(Flux3[i]/Flux4[i])

sResult1 = []
sResult2 = []
sResult3 = []
for i in range(0,numStars):
    sResult1.append(Sigma1[i]/Sigma4[i])
    sResult2.append(Sigma2[i]/Sigma4[i])
    sResult3.append(Sigma3[i]/Sigma4[i])

print(Flux4)
plt.plot(Result1, 'o')
plt.plot(Result2, 'o')
plt.plot(Result3, 'o')
plt.legend(['Hann', 'Gaussian', 'Exact'])
plt.title('Ratio of Method/True flux')
plt.show()

plt.scatter(Result1, Flux4)
plt.scatter(Result2, Flux4)
plt.scatter(Result3, Flux4)
plt.legend(['Hann', 'Gaussian', 'Exact'])
plt.title("True Flux vs Ratio of methods")
plt.show()


plt.plot(sResult1, 'o')
plt.plot(sResult2, 'o')
plt.plot(sResult3, 'o')
plt.legend(['Hann', 'Gaussian', 'Exact'])
plt.title('Method/True Value - Sigma')
plt.show()

plt.scatter(sResult1, Sigma4)
plt.scatter(sResult2, Sigma4)
plt.scatter(sResult3, Sigma4)
plt.legend(['Hann', 'Gaussian', 'Exact'])
plt.title("True Sigma vs Ratio of methods")
plt.show()





