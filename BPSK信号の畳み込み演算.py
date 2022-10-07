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
plt.rcParams["font.size"] = 14

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
Nbits = 25             # number of bits
f0 = 30*10**9          # carrier Frequency
fg = 1*10**9           # 情報信号の繰り返し周波数
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
#inputBits = np.random.randn(Nbits,1) > 0
# 生成するbitパターンを固定化する
inputBits = np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0]) > 0
inputBits = inputBits.reshape((25,1))

#はじめから1bitずらした信号を生成することで遅延を再現
inputBits2 = np.insert(inputBits, 0, 0, axis=0)
inputBits2 = np.delete(inputBits2, 25, 0)
print(inputBits.shape, inputBits2.shape)

#Digital-to-Analog Conversion
inputSignal = (np.tile(inputBits*2-1,(1,Ns))).ravel()
inputSignal2 = (np.tile(inputBits2*2-1,(1,Ns))).ravel()

print(inputSignal.shape, inputSignal2.shape)

#Multiplicator / mixxer
AM_signal = inputSignal*( carrier1)# + intermodulation1+ intermodulation2)
AM_signal2 = inputSignal2*( carrier1)


#----- 畳み込み演算 -----#
CONV = sig.convolve(AM_signal, AM_signal2, mode='same', method='fft')
CONV_normalized = CONV / np.max(CONV)
#yh = np.abs(sig.hilbert(AM_signal))

#----- 乗算（multiplication ） -----#
Multiple = AM_signal * AM_signal2
Multiple_normalized = Multiple / np.max(Multiple)


#---------- Plot of amplitude modulated signal ------------#
fig, axis = plt.subplots(5, 1)
fig.suptitle('BPSK Modulation')

axis[0].plot(t, inputSignal, color='C1')
axis[0].set_title('NRZ signal') # (Source Code/ Block Diagram: "inputSignal")
axis[0].set_xlabel('Time [s]')
axis[0].set_xlim(0,timeDomainVisibleLimit)
axis[0].set_ylabel('Amplitude [V]')
axis[0].grid(linestyle='dotted')

'''
axis[1].plot(t, inputSignal2, color='C2')
axis[1].set_title('NRZ signal (1ns delay)') # (Source Code/ Block Diagram: "carrier1")
axis[1].set_xlabel('Time [s]')
axis[1].set_xlim(0,timeDomainVisibleLimit)
axis[1].set_ylabel('Amplitude [V]')
axis[1].grid(linestyle='dotted')
'''


axis[1].plot(t,AM_signal, color='C3')
axis[1].set_title('BPSK signal') # (Source Code/ Block Diagram: "BPSK_signal")
axis[1].set_xlabel('Time [s]')
axis[1].set_xlim(0,timeDomainVisibleLimit)
axis[1].set_ylabel('Amplitude [V]')
axis[1].grid(linestyle='dotted')

axis[2].plot(t, AM_signal2, color='C2')
axis[2].set_title('BPSK signal (1ns delay)') # (Source Code/ Block Diagram: "BPSK_signal")
axis[2].set_xlabel('Time [s]')
axis[2].set_xlim(0,timeDomainVisibleLimit)
axis[2].set_ylabel('Amplitude [V]')
axis[2].grid(linestyle='dotted')

axis[3].plot(t, CONV_normalized, color='C4')
axis[3].set_xlabel('Time [s]')
axis[3].set_xlim(0,timeDomainVisibleLimit)
axis[3].set_ylabel('Intensity')
axis[3].set_title('Convolution')
axis[3].grid(linestyle='dotted')

axis[4].plot(t, Multiple_normalized, color='C4')
axis[4].set_xlabel('Time [s]')
axis[4].set_xlim(0,timeDomainVisibleLimit)
axis[4].set_ylabel('Intensity')
axis[4].set_title('Multiplication')
axis[4].grid(linestyle='dotted')

plt.subplots_adjust(hspace=0.85)

plt.show()