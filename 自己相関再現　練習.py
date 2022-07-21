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

from scipy.fftpack import fft
from numpy import e, pi, real
from numpy import sqrt
from numpy import sin
from numpy import cos
from numpy import zeros
from numpy import r_
from scipy.io.wavfile import read as wavread
from sympy import E, I

# 強度変調信号の生成
plt.rcParams["font.size"] = 18

# Used for symbol creation. Returns a decimal number from a 1 bit input
def GetBpskSymbol(bit1:bool):
    if(~bit1):
        return 0
    elif(bit1):
        return 1
    else:
        return -1

    # Maps a given symbol to a complex signal. Optionally, noise and phase offset can be added.
def BpskSymbolMapper(symbols:int,amplitude,noise1=0, noise2=0,  phaseOffset = 0):
    if(symbols == 0):
        return amplitude*(cos(np.deg2rad(0) + phaseOffset)) + (noise1 + 1j*noise2)
    elif(symbols == 1):
        return amplitude*(cos(np.deg2rad(180) + phaseOffset)) + (noise1 + 1j*noise2)
    else:
        return complex(0)

#-------------------------------------#
#---------- Configuration ------------#
#-------------------------------------#
fs = 1 * 10 ** 12           # sampling rate
baud = 1 * 10 ** 9          # symbol rate = bps?
Nbits = 25                  # number of bits
f0 = 30 * 10 ** 9          # carrier Frequency
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

carrier1 = sin(2 * pi * f0 * t)

#----------------------------#
#---------- BPSK ------------#
#----------------------------#

# Modulator Input
inputBits = np.random.randn(Nbits,1) > 0

#Digital-to-Analog Conversion
inputSignal = (np.tile(inputBits,(1,Ns))).ravel()
dataSymbols = np.array([[GetBpskSymbol(inputBits[x])] for x in range(0,inputBits.size)])

#Multiplicator / mixxer
AM_signal = inputSignal*( carrier1)# + intermodulation1+ intermodulation2)

#---------- Preperation BPSK Constellation Diagram ------------#

amplitude = 1

#Generate noise. Two sources for uncorelated noise.
noiseStandardDeviation = 0.12
noise1 = np.random.normal(0,noiseStandardDeviation,dataSymbols.size)
noise2 = np.random.normal(0,noiseStandardDeviation,dataSymbols.size)

#---------- Plot of amplitude modulated signal ------------#
fig, axis = plt.subplots(3, 1)
fig.suptitle('Modulation')

axis[0].plot(t, inputSignal, color='C1')
axis[0].set_title('NRZ signal') # (Source Code/ Block Diagram: "inputSignal")
axis[0].set_xlabel('Time [s]')
axis[0].set_xlim(0,timeDomainVisibleLimit)
axis[0].set_ylabel('Amplitude [V]')
axis[0].grid(linestyle='dotted')

axis[1].plot(t, carrier1, color='C2')
axis[1].set_title('Carrier signal') # (Source Code/ Block Diagram: "carrier1")
axis[1].set_xlabel('Time [s]')
axis[1].set_xlim(0,timeDomainVisibleLimit)
axis[1].set_ylabel('Amplitude [V]')
axis[1].grid(linestyle='dotted')

axis[2].plot(t,AM_signal, color='C3')
axis[2].set_title('Amplitude modulated signal') # (Source Code/ Block Diagram: "BPSK_signal")
axis[2].set_xlabel('Time [s]')
axis[2].set_xlim(0,timeDomainVisibleLimit)
axis[2].set_ylabel('Amplitude [V]')
axis[2].grid(linestyle='dotted')

plt.subplots_adjust(hspace=0.6)


#----- 相互相関関数の練習 -----#
ACF = np.correlate(AM_signal, AM_signal, mode="same")
print(ACF)
# 正規化して最大値±１となるようにグラフを描いてみる

fig, ax = plt.subplots()
ax.plot(t, ACF, color='C1')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Intensity')
#fig.suptitle('BPSK Modulation', fontsize=18)

# ax.set_title('Magnitude Spectrum (Source Code/ Block Diagram: "BPSK_signal")')

plt.show()