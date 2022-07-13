# Import functions and libraries
import numpy as np
import dask.array as da
import matplotlib.pyplot as plt
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

plt.rcParams["font.size"] = 18

# シンボル生成、1bitの入力から十進数を返す
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
symbolsToShow = 25          # グラフに表示される範囲 = bit数
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
inputSignal = (np.tile(inputBits*2-1,(1,Ns))).ravel() #int32
dataSymbols = np.array([[GetBpskSymbol(inputBits[x])] for x in range(0,inputBits.size)]) #int32

#Multiplicator / mixxer
BPSK_signal = inputSignal*(carrier1) # + intermodulation1+ intermodulation2) float64
#BPSK_signal = BPSK_signal.astype(np.float16)

#---------- 信号の積分（光信号に変換？） ------------#
f1, t1 = sp.symbols('f1, t1')

fx = sp.exp(-sp.I * 2 * sp.pi * f1 * t1) #BPSK_signalを掛けるとエラーが発生し積分できない
integ = sp.integrate(fx, (t1, 0, 25*10**-9), conds = 'none')

print(integ)
print(integ.subs(f1,100))
integ = (integ * BPSK_signal)
print(integ)

print(integ.dtype)

#f1 =np.arange(-12500, 12500)
print(integ[250].subs(f1, 100))

p1 = sp.plot(integ, (f1, -12500, 12500))
