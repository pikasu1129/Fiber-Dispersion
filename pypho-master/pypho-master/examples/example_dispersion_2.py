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
from pypho_functions import *
import numpy as np
import copy
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# Define some parameter
T_FWHM = 50e-12       # FWHM-pulsewidth of the optical power
T_e = T_FWHM / 1.66510  # 1/e-pulsewidth of the electrical field
z = 29.756e3         # Fiber length
D = 16.8                # Dispersion coefficient
T_FWHM_norm_max = 0.75   # Maximum FWHM pulseiwdth of the 

# Define network elements
gp       = pypho_setup(nos = 4*16, sps = 1*256, symbolrate = 10e9)
bitsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'random') # Set pattern = "singlepulse"
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = T_FWHM*gp.symbolrate)
sig_1550 = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = 0)
SSMF     = pypho_fiber(glova = gp, l = z,  D = D,   S = 0, alpha = 0.2e-12, gamma = 0, phi_max = 10.0)
eye      = pypho_eye(glova = gp, polarisation = 'x,y')

# Simulation
bits = bitsrc()

esig = esigsrc(bitsequence = bits)
E_Tx = sig_1550(esig = esig)                                              

E = copy.deepcopy(E_Tx)

# Fiber trannsmission
E = SSMF(E = E) 
E[0]['E'][0] = E[0]['E'][0] /np.max(np.abs(E[0]['E'][0]))*1e-3
E_Tx[0]['E'][0] = E_Tx[0]['E'][0] /np.max(np.abs(E_Tx[0]['E'][0]))*1e-3
# Plot Input and Output signal
beta_2, beta_3 = DS2beta(D, 0, gp.lambda0)
L_D = (T_e)**2  / np.abs(beta_2)     # Pulse width T_e here of the electrical field
 
plt.figure(1)
eye(E = E_Tx[0]['E'], polarisation = 'x', style="r")
eye(E = E[0]['E'], polarisation = 'x', style="g")

red_patch = mpatches.Patch(color='red', label='Input signal')
green_patch = mpatches.Patch(color='green', label='Output signal')
plt.legend(handles=[red_patch, green_patch], loc=1); plt.show()

print ('Calculated output FWHM-pulse width : T_FWHM,1 = ' , T_FWHM * np.sqrt(1 + (z/L_D)**2)*1e12, ' ps')
T_FWHM_norm_0 = T_FWHM*1e12 / 100.0
print (  'Calculated max reach : ' , 100e-12**2 *(1/1.665095)**2 * 2.0*np.pi * 299792458 *  T_FWHM_norm_0*np.sqrt(  T_FWHM_norm_max**2 - T_FWHM_norm_0**2) / (16.8e-6 * 1550e-9**2 )  *1e-3 )