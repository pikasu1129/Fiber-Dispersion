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
from pypho_functions import *
import numpy as np
import copy
import matplotlib.pyplot as plt


# Define network elements
gp       = pypho_setup(nos = 32, sps = 1*128, symbolrate = 10e9)
bitsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'ones')
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.33)
sig      = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = np.pi/8)
SSMF     = pypho_fiber(glova = gp, l = 80e3,  D = 0.0,   S = 0.0, alpha = 0.2, gamma = 1.4, phi_max = .01)

# Simulation

# Define wavelength channel
bits    = bitsrc()
esig    = esigsrc(bitsequence = bits)
E       = sig(esig = esig)                                              

E_Tx = copy.deepcopy(E)


# Fiber transmission
       
E = SSMF(E = E)    

# Plot power and phase of both pol axis
plt.figure(1)
plt.subplot(2, 1, 1)
plt.plot(gp.timeax()*1.0e12, np.abs(E[0]['E'][0])**2, 'r', label='$E_x(0, t)$')
plt.plot(gp.timeax()*1.0e12, np.abs(E[0]['E'][1])**2, 'g', label='$E_y(0, t)$')
plt.ylabel('$10log |E_{x,y}|^2$'); plt.xlabel('Time [ps]');
plt.grid(True)
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(gp.timeax()*1.0e12, np.angle(E[0]['E'][0]), 'r', label='$E_x(0, t)$')
plt.plot(gp.timeax()*1.0e12, np.angle(E[0]['E'][1]), 'g', label='$E_y(0, t)$')
plt.ylabel('$ \phi_{x,y} $'); plt.xlabel('Time [ps]');
plt.grid(True)
plt.legend()
plt.show()

# Calculate phase shift
L_eff = (1-np.exp(-SSMF.l*SSMF.alpha))/(SSMF.alpha)
phi_x_max = - (np.max(np.abs(E_Tx[0]['E'][0])**2) + 2.0/3.0* np.max(np.abs(E_Tx[0]['E'][1])**2))* L_eff * SSMF.gamma *1e-3
phi_y_max = - (np.max(np.abs(E_Tx[0]['E'][1])**2) + 2.0/3.0* np.max(np.abs(E_Tx[0]['E'][0])**2))* L_eff * SSMF.gamma *1e-3

print ('L_eff = ', L_eff)
print ('phi_x_max = ', phi_x_max)
print ('phi_y_max = ', phi_y_max)