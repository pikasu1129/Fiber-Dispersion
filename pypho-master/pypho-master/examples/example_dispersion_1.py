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
from scipy.interpolate import UnivariateSpline

# Define network elements
gp       = pypho_setup(nos =2**4, sps = 256, symbolrate = 10.0e9)


symbolsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'debruijn')
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'rect' , fwhm = 0.85)
sig_1550 = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = 0)
SSMF     = pypho_fiber(glova = gp, l = 60.0e3,  D = 17.0,   S = 0, alpha = 0.2e-12, gamma = 1.4e-12, phi_max = 10.0)

# Simulation
bits = symbolsrc()
esig = esigsrc(bitsequence = bits)

E_Tx = sig_1550(esig = esig)                                                      



# Define your parameters here
T_0 = 25.0e-12
z = SSMF.l
D = 17.0
beta_2, beta_3 = DS2beta(17.0, 0, gp.lambda0)


# Create a single pulse with gaussian shape (not power!)
E_Tx[0]['E'][0] = E_Tx[0]['E'][0]*0 + np.exp(-(gp.timeax()-gp.timeax()[-1]/2)**2 / (2.0*T_0**2) )

E = copy.deepcopy(E_Tx)


# Fiber transmission
E = SSMF(E = E, D = D, l = z) 
#sys.exit()


# Plot Input and Output signal 
plt.figure(1)
plt.plot(gp.timeax()*1.0e12, np.abs(E_Tx[0]['E'][0]), 'r', label='$E(0, t)$')
plt.plot(gp.timeax()*1.0e12, np.abs(E[0]['E'][0]), 'g', label= '$E(z='+ str(SSMF.l/1e3)+ 'km, t)')

# Get FWHM of the input signal E_Tx
spline = UnivariateSpline(gp.timeax()*1.0e12, np.abs(E_Tx[0]['E'][0])-1*np.max(np.abs(E_Tx[0]['E'][0]))/2, s=0)
r1, r2 = spline.roots() # find the roots
plt.annotate(s='', xy=(r1,np.max(np.abs(E_Tx[0]['E'][0]))/2), xytext=(r2,np.max(np.abs(E_Tx[0]['E'][0]))/2), arrowprops=dict(arrowstyle='<->'))
plt.text(r1+(r2-r1)/2.0, 0.01 +np.max(np.abs(E_Tx[0]['E'][0]))/2, '$T_{FWHM,0}$ = ' + str(np.round(r2-r1,2)) + ' ps', fontsize=12, horizontalalignment='center')
T_FWHM_0 = (r2-r1) * 1e-12
T_0_plot = T_FWHM_0 / 2.35482


# Get FWHM of the output signal E
spline = UnivariateSpline(gp.timeax()*1.0e12, np.abs(E[0]['E'][0])-1*np.max(np.abs(E[0]['E'][0]))/2, s=0)
r1, r2 = spline.roots() # find the roots
plt.annotate(s='', xy=(r1,np.max(np.abs(E[0]['E'][0]))/2), xytext=(r2,np.max(np.abs(E[0]['E'][0]))/2), arrowprops=dict(arrowstyle='<->'))
plt.text(r1+(r2-r1)/2.0, 0.01 + np.max(np.abs(E[0]['E'][0]))/2, '$T_{FWHM,1}$ = ' + str(np.round(r2-r1,2)) + ' ps', fontsize=12, horizontalalignment='center')
T_FWHM_1 = (r2-r1) * 1e-12
plt.ylabel('$|E|$ a.u.'); plt.xlabel('Time $t$ [ps]'); legend = plt.legend(loc='upper right')
plt.show()
L_D = (T_0_plot)**2  / np.abs(beta_2)
# Print the results
print ('Input signal 1/e-pulse width by definition : T_0 = ' , T_0*1e12, ' ps')
print ('Input signal 1/e-pulse width from plot : T_0 = ' , T_0_plot*1e12, ' ps')
print ('Input signal FWHM-pulse width from plot : T_FWHM,0 = ' , T_FWHM_0*1e12, ' ps')
print ('Output signal FWHM-pulse width from plot : T_FWHM,1 = ' , T_FWHM_1*1e12, ' ps')
print( 'Calculated output FWHM-pulse width : T_FWHM,1 = ' , T_FWHM_0 * np.sqrt(1 + (z/L_D)**2)*1e12, ' ps')