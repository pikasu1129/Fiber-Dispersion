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
plt.rcParams["font.size"] = 17

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
'''
# Modulator Input
#inputBits = np.random.randn(Nbits,1) > 0
# 生成するbitパターンを固定化する
inputBits = np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0]) > 0
inputBits = inputBits.reshape((25,1))

#はじめから1bitずらした信号を生成することで遅延を
inputBits2 = np.insert(inputBits, 0, 0, axis=0)
inputBits2 = np.delete(inputBits2, 25, 0)
print(inputBits.shape, inputBits2.shape)

#Digital-to-Analog Conversion
inputSignal = (np.tile(inputBits,(1,Ns))).ravel()
inputSignal2 = (np.tile(inputBits2,(1,Ns))).ravel()

print(inputSignal.shape, inputSignal2.shape)
'''

#Multiplicator / mixxer
AM_signal = 0.5*(1 + 1*signal)*cos(2*pi*f0*t)
# 手計算で展開
AM_signal2 = 0.5*(cos(2*pi*f0*t) + 0.5*(cos(2*pi*(fg+f0)*t) + cos(2*pi*(fg-f0)*t)))

#AM_signal2 = inputSignal2*( carrier1)


#----- 相互相関関数or畳み込みの練習 -----#
CONV = np.convolve(AM_signal, AM_signal, mode="same") # fullにして時間軸をlen(ACF)にしてみる
CONV_normalized = CONV / np.max(CONV)
# 1ns遅延を与えてみる

#----- 乗算（multiplication ） -----#
Multiple = AM_signal * AM_signal
Multiple_normalized = Multiple / np.max(Multiple)

'''
ACF = np.convolve(AM_signal, AM_signal2, mode='same')

ACF_normalized = ACF / np.amax(ACF)
'''

#yh = np.abs(sig.hilbert(AM_signal))
# 正規化して最大値±１となるようにグラフを描いてみる

#---------- Plot of amplitude modulated signal ------------#
fig, axis = plt.subplots(4, 1)
fig.suptitle('Modulation')

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

axis[3].plot(t, CONV_normalized, color='C4')
axis[3].set_xlabel('Time [s]')
axis[3].set_xlim(0,timeDomainVisibleLimit)
axis[3].set_ylabel('Intensity')
axis[3].set_title('Convolution')
axis[3].grid(linestyle='dotted')


'''
axis[3].plot(t, AM_signal2, color='C2')
axis[3].set_title('Amplitude modulated signal (1ns delay)') # (Source Code/ Block Diagram: "BPSK_signal")
axis[3].set_xlabel('Time [s]')
axis[3].set_xlim(0,timeDomainVisibleLimit)
axis[3].set_ylabel('Amplitude [V]')
axis[3].grid(linestyle='dotted')
'''

#---------- Plot of convvolution function ------------#
fig = plt.figure()
ax = fig.add_subplot(211)
ax.plot(t, CONV, color='C4')
#ax.set_xlabel('Time [s]')
ax.set_ylabel('Intensity')
ax.set_title('Convolution')
#fig.suptitle('Auto  Correlation Function')

ax2 = fig.add_subplot(212)
ax2.plot(t, Multiple_normalized, color='C1')
ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Intensity')
ax2.set_title('Multiplication')

#---------- Plot of Modulated Signal and Spectrum ------------#
fig = plt.figure(constrained_layout=True)
gs = gridspec.GridSpec(4, 1, figure=fig)
fig.suptitle('BPSK Modulation', fontsize=18)
ax = fig.add_subplot(gs[0, :])

ax1 = fig.add_subplot(gs[0])
ax1.set_title('Magnitude Spectrum (Source Code/ Block Diagram: "AM_signal")')
ax1.magnitude_spectrum(AM_signal2, Fs=fs, color='C1')
ax1.set_xlim(lowerLimit,upperLimit)
ax1.grid(linestyle='dotted')

ax2 = fig.add_subplot(gs[1])
ax2.set_title('Log. Magnitude Spectrum (Source Code/ Block Diagram: "AM_signal")')
ax2.magnitude_spectrum(AM_signal2, Fs=fs, scale='dB', color='C1')
ax2.set_xlim(lowerLimit,upperLimit)
ax2.grid(linestyle='dotted')


ax3 = fig.add_subplot(gs[2])
ax3.set_title('Power Spectrum Density (PSD) (Source Code/ Block Diagram: "AM_signal")')
ax3.psd(AM_signal2,NFFT=len(t),Fs=fs)
ax3.set_xlim(lowerLimit,upperLimit)
ax3.grid(linestyle='dotted')


ax4 = fig.add_subplot(gs[3])
ax4.set_title('Carrier")')
ax4.magnitude_spectrum(carrier1, Fs=fs, color='C1')
ax4.set_xlim(lowerLimit,upperLimit)
ax4.grid(linestyle='dotted')


plt.subplots_adjust(hspace=0.85)

plt.show()