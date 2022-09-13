##!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Import functions and libraries
import sys
sys.path.append('../')
from pypho_setup import pypho_setup
from pypho_symbols import pypho_symbols
from pypho_signalsrc import pypho_signalsrc
from pypho_lasmod import pypho_lasmod
from pypho_cfiber import pypho_fiber
from pypho_cwlaser import pypho_cwlaser
from pypho_optfi import pypho_optfi
from pypho_functions import *
import numpy as np
import copy
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt



# Define network elements
gp       = pypho_setup(nos = 8, sps = 256, symbolrate = 10e9)
bitsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'singlepulse') # Set pattern = "singlepulse"
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.33)
sig_f0   = pypho_cwlaser(glova = gp, power = -20, Df = 0,  teta =1*np.pi/4.0)
sig_f1   = pypho_lasmod(glova = gp, power = -10, Df = 500, teta = 0*np.pi/1.0)
filter_f0 = pypho_optfi(glova = gp, Df = 0, B = 150)
SSMF     = pypho_fiber(glova = gp, l = 10e3,  D = 0.0,   S = 0.0, alpha = 0.2, gamma = 1.4, phi_max = .001)

# Simulation

# Define wavelength channels

E_f0    = sig_f0()                                              

bits    = bitsrc()
esig    = esigsrc(bitsequence = bits)
E_f1    = sig_f1(esig = esig)  
E       = copy.deepcopy(E_f1)

E[0]['E'][0] = E_f0[0]['E'][0] + E_f1[0]['E'][0]  # Multiplex all signals X-Pol
E[0]['E'][1] = E_f0[0]['E'][1] + E_f1[0]['E'][1]  # Multiplex all signals Y-Pol 


E_Tx = copy.deepcopy(E)
# Plot Spectrum Input signals
plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot((gp.freqax()-gp.f0)*1e-9, np.log10(abs(fftshift(fft(E[0]['E'][0]   )**2))), 'r', label='X-Pol')
plt.plot((gp.freqax()-gp.f0)*1e-9, np.log10(abs(fftshift(fft(E[0]['E'][1]   )**2))), 'g', label='Y-Pol')
plt.title("Input spectrum", loc='left')
plt.ylabel('Spec. density');
plt.grid(True)
plt.legend()
# Fiber transmission
       
E = SSMF(E = E)    

# Plot Spectrum Output signals
plt.subplot(2, 1, 2)
plt.plot((gp.freqax()-gp.f0)*1e-9, np.log10(abs(fftshift(fft(E[0]['E'][0]   )**2))), 'r', label='X-Pol')
plt.plot((gp.freqax()-gp.f0)*1e-9, np.log10(abs(fftshift(fft(E[0]['E'][1]   )**2))), 'g', label='Y-Pol')
plt.title("Output spectrum", loc='left')
plt.ylabel('Spec. density'); plt.xlabel('Frequency deviation [GHz]');
plt.grid(True)


# Calculate phase shift
L_eff = (1-np.exp(-SSMF.l*SSMF.alpha))/(SSMF.alpha)
phi_XPM_x_max = - 2*(np.max(np.abs(E_f1[0]['E'][0])**2) + 2.0/3.0* np.max(np.abs(E_f1[0]['E'][1])**2))* L_eff * SSMF.gamma *1e-3
phi_XPM_y_max = - 2*(np.max(np.abs(E_f1[0]['E'][1])**2) + 2.0/3.0* np.max(np.abs(E_f1[0]['E'][0])**2))* L_eff * SSMF.gamma *1e-3

phi_SPM_x_max = - (np.max(np.abs(E_f0[0]['E'][0])**2) + 2.0/3.0* np.max(np.abs(E_f0[0]['E'][1])**2))* L_eff * SSMF.gamma *1e-3
phi_SPM_y_max = - (np.max(np.abs(E_f0[0]['E'][1])**2) + 2.0/3.0* np.max(np.abs(E_f0[0]['E'][0])**2))* L_eff * SSMF.gamma *1e-3

print ('L_eff = ', L_eff)
print ('phi_XPM_x_max = ', phi_XPM_x_max)
print ('phi_XPM_y_max = ', phi_XPM_y_max)
print ('phi_SPM_x_max = ', phi_SPM_x_max)
print ('phi_SPM_y_max = ', phi_SPM_y_max)
print ('phi_x_max = phi_XPM_x_max + phi_SPM_x_max', phi_XPM_x_max + phi_SPM_x_max)
print ('phi_y_max = phi_XPM_y_max + phi_SPM_y_max', phi_XPM_y_max + phi_SPM_y_max)

# Plot power of both pol axis as function of transmission distance

E = filter_f0(E)

plt.figure(2)

plt.subplot(2, 1, 1)
plt.ylabel('$|E|^2$'); plt.xlabel('Transmission distance [m]');
plt.plot(gp.timeax()*1.0e12, np.abs(E_f1[0]['E'][0])**2, 'r', label='$Pulse: E_x $')
plt.plot(gp.timeax()*1.0e12, np.abs(E_f1[0]['E'][1])**2, 'g', label='$Pulse: E_y Pulse$')
plt.plot(gp.timeax()*1.0e12, np.abs(E_f0[0]['E'][0])**2, 'r:', label='$cw@1550nm: E_x $')
plt.plot(gp.timeax()*1.0e12, np.abs(E_f0[0]['E'][1])**2, 'g:', label='$cw@1550nm: E_y$')
legend = plt.legend(loc='upper right')

plt.grid(True)
plt.subplot(2, 1, 2)
plt.plot(gp.timeax()*1.0e12, np.angle(E[0]['E'][0]), 'r', label='$cw@1550nm: E_x$')
plt.plot(gp.timeax()*1.0e12, np.angle(E[0]['E'][1]), 'g', label='$cw@1550nm: E_y$')
plt.ylabel('$phase(E_x)$');
plt.grid(True)
legend = plt.legend(loc='upper right')
plt.show()