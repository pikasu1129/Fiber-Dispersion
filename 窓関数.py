import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

from scipy import signal
from numpy import r_

#-------------------------------------#
#---------- Configuration ------------#
#-------------------------------------#
fs = 1 * 10 ** 12           # sampling rate
baud = 1 * 10 ** 9          # symbol rate = bps?
Nbits = 25                  # number of bits
f0 = 30 * 10 ** 9           # carrier Frequency
Ns = int(fs/baud)           # number of Samples per Symbol
N = Nbits * Ns              # Total Number of Samples
t = r_[0.0:N]/fs            # time points
f = r_[0:N/2.0]/N*fs        # Frequency Points

w1 = np.hanning(N)
w2 = np.hamming(N)
w3 = signal.blackman(N)


fig = plt.figure(figsize=(7.0, 5.0))
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 12

plt.plot(w1, "r", w2, "b", w3, "g")
plt.axis("tight")
plt.ylabel("Amplitude")
plt.xlabel("data number")
plt.legend(["hanning", "hamming", "blackman"])

plt.show()