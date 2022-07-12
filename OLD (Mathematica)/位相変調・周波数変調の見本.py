import scipy.signal as sg
import numpy as np
import matplotlib.pyplot as plt

sample_rate = 44100.0
nsamples = 320
F_1 = 440.0
F_2 = 5000
t = np.arange(nsamples) / sample_rate
vin = np.sin(2 * np.pi * F_1 * t)
vfm = np.sin(2 * np.pi * F_2 * t + 6.0 * -np.cos(2 * np.pi * F_1 * t))
vpm = np.sin(2 * np.pi * F_2 * t + 6.0 * -np.sin(2 * np.pi * F_1 * t))
fig = plt.figure(1)
ax = fig.add_subplot(311)
ax.plot(vin[1: 300])
ax = fig.add_subplot(312)
ax.plot(vfm[1: 300])
ax = fig.add_subplot(313)
ax.plot(vpm[1: 300])
fig.set_tight_layout(True)
plt.savefig("fm2.png")
plt.show()