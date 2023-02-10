import numpy as np
import matplotlib.pyplot as plt

#antiquated

file = 'Studentship/WD0830_Results/RMSResults.txt'
noise = []
time = []
sqrtTime = []
with open(file, 'r') as f:
    lines = f.readlines()
    for line in lines:
        var = line.strip('\n').split(',')
        time.append(float(var[1]))
        noise.append(float(var[0]))
for item in time:
    sqrtTime.append(np.sqrt(item))

plt.scatter(time, noise)
plt.title('Noise vs Time, WD0830 B Filter')
plt.xlabel('Time (seconds)')
plt.ylabel("Noise (pixels)")
plt.show()
plt.savefig("Studentship/WD0830_Results/NoiseVsTime.png")

plt.scatter(sqrtTime, noise)
plt.title('Noise vs sqrt(Time), WD0830 B Filter')
plt.xlabel('Time (seconds)')
plt.ylabel("Noise (pixels)")
plt.show()
plt.savefig("Studentship/WD0830_Results/NoiseVsSqrtTime.png")