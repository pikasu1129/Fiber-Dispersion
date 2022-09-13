##!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Import functions and libraries
import sys
sys.path.append('../')
from pypho_setup import pypho_setup
from pypho_cfiber import pypho_fiber
from pypho_cwlaser import pypho_cwlaser
from pypho_optfi import pypho_optfi
from pypho_functions import *
import numpy as np
import copy
import matplotlib.pyplot as plt
import time
#  Define some paramaters

delta_f = 25       # Frequency difference of the two cw-signals
L = 80e3           # Fiber length in m
l = 0.2e3           # Step fiber length in m
k = 7               # Defines how man channels left and right will be analyzed
chut = np.arange(-k*delta_f, k*delta_f+delta_f, delta_f)        # Channels under test


P_FWM_x = np.zeros((chut.size, int(L/l)))
P_FWM_y = np.zeros((chut.size, int(L/l)))

# Define network elements
gp       = pypho_setup(nos = 2**5, sps = 128, symbolrate = 10e9)
sig_f0  = pypho_cwlaser(glova = gp, power = 3, Df = 0, teta = 1*np.pi/4.0)
sig_f   = pypho_cwlaser(glova = gp, power = 3, Df = delta_f, teta = 1*np.pi/2.0)
filter_df = pypho_optfi(glova = gp, Df = 0, B = 10)
SSMF     = pypho_fiber(glova = gp, l = l,  D = 5.0,   S = 0.0, alpha = 0.2e-12, gamma = 1.4, phi_max = .1)

# Simulation

# Define wavelength channel
E_f0  = sig_f0()                                              
E_fm1 = sig_f(Df = +delta_f, teta = np.pi/2)
E_fp1 = sig_f(Df = -delta_f, teta = np.pi/2)
E_fm2 = sig_f(Df = -2*delta_f, teta = np.pi/2)
E_fp2 = sig_f(Df = +2*delta_f, teta = np.pi/2)
E = copy.deepcopy(E_fm1)

E[0]['E'][0] = E_f0[0]['E'][0] + E_fp1[0]['E'][0] + E_fm1[0]['E'][0] + E_fp2[0]['E'][0] + E_fm2[0]['E'][0] # Multiplex all signals X-Pol
E[0]['E'][1] = E_f0[0]['E'][1] + E_fp1[0]['E'][1] + E_fm1[0]['E'][1] + E_fp2[0]['E'][1] + E_fm2[0]['E'][1] # Multiplex all signals Y-Pol 


E_Tx = copy.deepcopy(E)
# Plot Spectrum Input signals
plt.figure(1)
plt.subplot(2, 1, 1)
plt.semilogy((gp.freqax()-gp.f0)*1e-9, abs(fftshift(fft(E[0]['E'][0]   )))**2, 'r')
plt.semilogy((gp.freqax()-gp.f0)*1e-9, abs(fftshift(fft(E[0]['E'][1]   )))**2, 'g:')
plt.title("Input spectrum", loc='left')
plt.ylabel('Spec. density');
plt.grid(True)
# Fiber transmission
c1 = 0
t=time.time()
for z in np.arange(0,L,l):
    print(z)

    c2 = 0
    for ch in chut:         # Loop through wavelengths and save mean power value of both pol axis
        E_in = copy.deepcopy(E)

        E_filtered = filter_df(E_in, Df = ch)

        P_FWM_x[c2][c1] = np.mean( np.abs(E_filtered[0]['E'][0])**2 )
        P_FWM_y[c2][c1] = np.mean( np.abs(E_filtered[0]['E'][1])**2 )
        c2 = c2 + 1

    c1 = c1 + 1        
    E = SSMF(E = E)    
print('Time='+str(time.time()-t))
# Plot Spectrum Output signals
plt.subplot(2, 1, 2)
plt.semilogy((gp.freqax()-gp.f0)*1e-9, abs(fftshift(fft(E[0]['E'][0]   )))**2, 'r')
plt.semilogy((gp.freqax()-gp.f0)*1e-9, abs(fftshift(fft(E[0]['E'][1]   )))**2, 'g:')
plt.title("Output spectrum", loc='left')
plt.ylabel('Spec. density'); plt.xlabel('Frequency deviation [GHz]');
plt.grid(True)

# Plot power of both pol axis as function of transmission distance
plt.figure(2)
for i in np.arange(0,chut.size):
    plt.subplot(2, 1, 1)
    plt.semilogy(np.arange(0,L,l), P_FWM_x[i]*1e3, label=str(chut[i]) + ' GHz')
    plt.subplot(2, 1, 2)
    plt.semilogy(np.arange(0,L,l), P_FWM_y[i]*1e3)

plt.subplot(2, 1, 2)
plt.ylabel('$10log |E_y|^2$'); plt.xlabel('Transmission distance [m]');
plt.grid(True)
plt.subplot(2, 1, 1)
plt.ylabel('$10log |E_y|^2$');
plt.grid(True)
plt.legend()
plt.show()