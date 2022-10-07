# Import functions and libraries
import scipy.signal as sig
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from IPython.core.display import display
import matplotlib.gridspec as gridspec
import scipy as scy
import sympy as sp
import threading,time
import multiprocessing
import sys
import bitarray
import cmath

from scipy.fftpack import fft, ifft
from numpy import e, ndarray, pi, real
from numpy import sqrt
from numpy import sin
from numpy import cos
from numpy import zeros
from numpy import r_
from scipy.io.wavfile import read as wavread
from sympy import E, I, NDimArray

# 強度変調信号の生成
plt.rcParams["font.size"] = 10

plt.subplots_adjust(hspace=0.75)


# Used for symbol creation. Returns a decimal number from a 1 bit input
def GetBpskSymbol(bit1:bool):
    if(~bit1):
        return 0
    elif(bit1):
        return 1
    else:
        return -1

#-------------------------------------#
#---------- Configuration ------------#
#-------------------------------------#
fs = 1*10**12          # sampling rate
baud = 1*10**9         # symbol rate = bps?
Nbits = 2**9-1                  # number of bits
f0 = 30*10**9          # carrier Frequency
fg = 1*10**9          # 情報信号の繰り返し周波数
Ns = int(fs/baud)           # number of Samples per Symbol
N = Nbits * Ns              # Total Number of Samples
t = r_[0.0:N]            # time points float64
f = r_[0:N/2.0]       # Frequency Points float64
t = t.astype(np.float32)
f = f.astype(np.float32)
t = t/fs
f = f/N*fs

# Limit for representation of time domain signals for better visibility.
symbolsToShow = 25
timeDomainVisibleLimit = np.minimum(Nbits/baud,symbolsToShow/baud)

# Limit for representation of frequency domain signals for better visibility.
sideLobesToShow = 9
sideLobeWidthSpectrum = baud
lowerLimit = np.maximum(0,f0-sideLobeWidthSpectrum*(1 + sideLobesToShow))
upperLimit = f0 + sideLobeWidthSpectrum*(1 + sideLobesToShow)

carrier1 = cos(2 * pi * f0 * t)
signal = cos(2 * pi * fg * t)

#----------------------------#
#---------- AM   ------------#
#----------------------------#

#Multiplicator / mixxer
AM_signal = 0.5*(1 + 1*signal)*cos(2*pi*f0*t)
# 手計算で展開
AM_signal2 = 0.5*(cos(2*pi*f0*t) + 0.5*(cos(2*pi*(fg+f0)*t) + cos(2*pi*(fg-f0)*t)))
AM_signal3 = 0.5*(cos(2*pi*f0*t) + 0.5*(cos(2*pi*(fg+f0)*t) - cos(2*pi*(fg-f0)*t)))
#AM_signal2 = inputSignal2*( carrier1)


# 正規化して最大値±１となるようにグラフを描いてみる

#---------- Plot of amplitude modulated signal ------------#
fig, axis = plt.subplots(4, 1)
fig.suptitle('Modulation')
fig.subplots_adjust(hspace=0.75)

axis[0].plot(t, signal, color='C1')
axis[0].set_title('signal') # (Source Code/ Block Diagram: "inputSignal")
axis[0].set_xlabel('Time [s]')
axis[0].set_xlim(0,timeDomainVisibleLimit)
axis[0].set_ylabel('Amplitude [V]')
axis[0].grid(linestyle='dotted')

axis[1].plot(t, carrier1, color='C2')
axis[1].set_title('carrier signal (1ns delay)') # (Source Code/ Block Diagram: "carrier1")
axis[1].set_xlabel('Time [s]')
axis[1].set_xlim(0,timeDomainVisibleLimit)
axis[1].set_ylabel('Amplitude [V]')
axis[1].grid(linestyle='dotted')

axis[2].plot(t,AM_signal2, color='C3')
axis[2].set_title('Amplitude modulated signal') # (Source Code/ Block Diagram: "BPSK_signal")
axis[2].set_xlabel('Time [s]')
axis[2].set_xlim(0,timeDomainVisibleLimit)
axis[2].set_ylabel('Amplitude [V]')
axis[2].grid(linestyle='dotted')


axis[3].plot(t, AM_signal3, color='C3')
axis[3].set_xlabel('Time [s]')
axis[3].set_xlim(0,timeDomainVisibleLimit)
axis[3].set_ylabel('Intensity')
axis[3].set_title('LSB*-1')
axis[3].grid(linestyle='dotted')


#---------- Plot of Modulated Signal and Spectrum ------------#
fig = plt.figure(constrained_layout=True)
gs = gridspec.GridSpec(4, 1, figure=fig)
fig.suptitle('BPSK Modulation', fontsize=18)
ax = fig.add_subplot(gs[0, :])

ax1 = fig.add_subplot(gs[1])
ax1.set_title('Magnitude Spectrum (Source Code/ Block Diagram: "AM_signal2")')
ax1.magnitude_spectrum(AM_signal2, Fs=fs, color='C1')
ax1.set_xlim(lowerLimit,upperLimit)
ax1.grid(linestyle='dotted')

ax2 = fig.add_subplot(gs[2])
ax2.set_title('Log. Magnitude Spectrum (Source Code/ Block Diagram: "AM_signal2")')
ax2.magnitude_spectrum(AM_signal2, Fs=fs, scale='dB', color='C1')
ax2.set_xlim(lowerLimit,upperLimit)
ax2.grid(linestyle='dotted')


ax3 = fig.add_subplot(gs[0])
ax3.set_title('Magnitude Spectrum (Source Code/ Block Diagram: "AM_signal3")')
ax3.magnitude_spectrum(AM_signal3, Fs=fs, color='C1')
ax3.set_xlim(lowerLimit,upperLimit)
ax3.grid(linestyle='dotted')


ax4 = fig.add_subplot(gs[3])
ax4.set_title('Log. Magnitude Spectrum (Source Code/ Block Diagram: "AM_signal3")')
ax4.magnitude_spectrum(AM_signal3, Fs=fs, scale='dB', color='C1')
ax4.set_xlim(lowerLimit,upperLimit)
ax4.grid(linestyle='dotted')

plt.show()