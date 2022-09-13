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
from pypho_eye import pypho_eye
from pypho_optfi import pypho_optfi
from pypho_functions import *
import numpy as np
import copy
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# Define network elements
gp       = pypho_setup(nos = 4*16, sps = 1*256, symbolrate = 40e9)
bitsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'random') # Set pattern = "singlepulse"
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.33)
sig_1550 = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = 0)
sig_1546 = pypho_lasmod(glova = gp, power = 0, Df = +500, teta = 0)
sig_1554 = pypho_lasmod(glova = gp, power = 0, Df = -500, teta = 0)
filter_1550 = pypho_optfi(glova = gp, Df = 0, B = 100)
filter_1546 = pypho_optfi(glova = gp, Df = +500, B = 100)
filter_1554 = pypho_optfi(glova = gp, Df = -500, B = 100)
SSMF     = pypho_fiber(glova = gp, l = 100e3,  D = 16.8,   S = 0.058, alpha = 0.2e-12, gamma = 0, phi_max = 10.0)
DCF      = pypho_fiber(glova = gp, l = 16.8e3,  D = -100.0,   S = -0.058/16.8 * 100.0, alpha = 0.2e-12, gamma = 0, phi_max = 10.0)
eye      = pypho_eye(glova = gp, polarisation = 'x,y')

# Simulation
bits = bitsrc()

esig = esigsrc(bitsequence = bits)
E_1550 = sig_1550(esig = esig)                                              
E_1546 = sig_1546(esig = esig)
E_1554 = sig_1554(esig = esig)


# Plot Input signals
plt.figure(1)
plt.subplot(3, 1, 1)
plt.title("$\lambda$=1546nm", loc='left')
eye(E = E_1546[0]['E'], polarisation = 'x', style="0.5")

plt.subplot(3, 1, 2)
plt.title("$\lambda$=1550nm", loc='left')
eye(E = E_1550[0]['E'], polarisation = 'x', style="0.5")

plt.subplot(3, 1, 3)
plt.title("$\lambda$=1554nm", loc='left')
eye(E = E_1554[0]['E'], polarisation = 'x', style="0.5")

E_Tx = copy.deepcopy(E_1550)
E_Tx[0]['E'][0] = E_1550[0]['E'][0] + E_1546[0]['E'][0] + E_1554[0]['E'][0] # Multiplex all signals 
E = copy.deepcopy(E_Tx)

# Fiber transmission
E = SSMF(E = E,l = 100e3,  D = 16.8, S = 0.058) 
E =  DCF(E = E,l = 16.8e3, D = -100, S = -0.058/16.8 * 100) # For full dispersion compensation set S = -0.058/16.8 * 100


# Filter out the channels
E_1550 = filter_1550(E = copy.deepcopy(E))
E_1546 = filter_1546(E = copy.deepcopy(E))
E_1554 = filter_1554(E = copy.deepcopy(E))

# Plot Output signals
#plt.figure(2)
plt.subplot(3, 1, 1)
eye(E = E_1546[0]['E'], polarisation = 'x', style="r")
red_patch = mpatches.Patch(color='0.5', label='Input signal')
green_patch = mpatches.Patch(color='red', label='Output signal')
plt.legend(handles=[red_patch, green_patch], loc=1)

plt.subplot(3, 1, 2)
eye(E = E_1550[0]['E'], polarisation = 'x', style="r")
plt.subplot(3, 1, 3)
eye(E = E_1554[0]['E'], polarisation = 'x', style="r")
plt.show()