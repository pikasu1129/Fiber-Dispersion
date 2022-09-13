##!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Import functions and libraries
import sys
sys.path.append('../')
from pypho_setup import pypho_setup
from pypho_symbols import pypho_symbols
from pypho_signalsrc import pypho_signalsrc
from pypho_lasmod import pypho_lasmod
from pypho_fiber import pypho_fiber
from pypho_fiber_birefringence import pypho_fiber_birefringence

import numpy as np
import copy

import matplotlib.pyplot as plt

# Define network elements
gp       = pypho_setup(nos = 16, sps = 1*128, symbolrate = 10e9)
bitsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'singlepulse') # Set pattern = "singlepulse"
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.33)
sig      = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = 1*np.pi/8.0)
SSMF     = pypho_fiber(glova = gp, l = 80e3,  D = 0.0,   S = 0.0, alpha = 1.0e-12, gamma = 1.4, phi_max = 0.9)

# Simulation

# Define wavelength channels
bits    = bitsrc()
esig    = esigsrc(bitsequence = bits)
E       = sig(esig = esig)  

E_Tx = copy.deepcopy(E)

# Create birefarray
fibres = []
fibres.append(pypho_fiber_birefringence(0, 0, 0.001))
#fibres.append(pypho_fiber_birefringence(40e3, np.pi/2, 0.001))

# Fiber transmission
E = SSMF(E = E, birefarray = fibres)    

# Plot power of both pol axis as function of transmission distance
plt.figure(2)
plt.ylabel('$|E|^2$'); plt.xlabel('Time [ps]');
plt.plot(gp.timeax()*1.0e12, np.abs(E_Tx[0]['E'][0])**2, 'r', alpha=0.4, label='$E_x$ input')
plt.plot(gp.timeax()*1.0e12, np.abs(E_Tx[0]['E'][1])**2, 'g', alpha=0.4, label='$E_y$ input')
plt.plot(gp.timeax()*1.0e12, np.abs(E[0]['E'][0])**2, 'r', label='$E_x$ output')
plt.plot(gp.timeax()*1.0e12, np.abs(E[0]['E'][1])**2, 'g', label='$E_y$ output')
plt.grid(True); legend = plt.legend(loc='upper right'); plt.show()