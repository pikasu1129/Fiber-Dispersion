##!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Import functions and libraries
# Relase 0.5
# arne.striegler@hm.edu

import sys
sys.path.append('../')
from pypho_setup import pypho_setup
from pypho_bits import pypho_bits
from pypho_signalsrc import pypho_signalsrc
from pypho_lasmod import pypho_lasmod
from pypho_fiber_mod import pypho_fiber
from pypho_functions import *
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# Define network elements
gp       = pypho_setup(nos = 16, sps = 256, symbolrate = 10e9)
bitsrc   = pypho_bits(glova = gp, nos = gp.nos, pattern = 'singlepulse')
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.85)
sig_1550 = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = 0)
SSMF     = pypho_fiber(glova = gp, l = 80.0e3,  D = 17,   S = 0, alpha = 0.2e-12, gamma = 0, phi_max=10.0)

# Simulation
bits = bitsrc()
esig = esigsrc(bitsequence = bits)
E_Tx = sig_1550(esig = esig)                                                      


# Define your parameters here
T_0 = 25.0e-12
z = 50.0e3
D = 17.0
beta_2, beta_3 = DS2beta(17.0, 0, gp.lambda0)


# Create a single pulse with gaussian shape (not power!)
E_Tx[0]['E'][0] = E_Tx[0]['E'][0]*0 + np.exp(-(gp.timeax-gp.timeax[-1]/2)**2 / (2.0*T_0**2) )
E = copy.deepcopy(E_Tx)
E[0]['E'][0] = E[0]['E'][0]*0 + np.sqrt(np.exp(-(gp.timeax-gp.timeax[-1]/2)**2 / (2.0*T_0**2) ) )

Eq = copy.deepcopy(E_Tx)
Eq[0]['E'][0] = Eq[0]['E'][0]*0 + (np.exp(-(gp.timeax-gp.timeax[-1]/2)**2 / (2.0*T_0**2) ) )**2


# Plot Input and Output signal 
plt.figure(1)
plt.plot(gp.timeax*1.0e12, np.abs(E_Tx[0]['E'][0]), 'r', label='Power: Gauss')
plt.plot(gp.timeax*1.0e12, np.abs(E[0]['E'][0]), 'g', label='Power: sqrt.Gauss')
plt.plot(gp.timeax*1.0e12, np.abs(Eq[0]['E'][0]), 'b', label='Power: Gauss**2')


# Get FWHM of the input signal E_Tx
spline = UnivariateSpline(gp.timeax*1.0e12, np.abs(E_Tx[0]['E'][0])-1*np.max(np.abs(E_Tx[0]['E'][0]))/2, s=0)
r1, r2 = spline.roots() # find the roots
plt.annotate(s='', xy=(r1,np.max(np.abs(E_Tx[0]['E'][0]))/2), xytext=(r2,np.max(np.abs(E_Tx[0]['E'][0]))/2), arrowprops=dict(arrowstyle='<->'))
plt.text(r1+(r2-r1)/2.0, 0.04 +np.max(np.abs(E_Tx[0]['E'][0]))/2, '$T_{FWHM_0}$ = ' + str(np.round(r2-r1,2)) + ' ps', fontsize=12, horizontalalignment='center')
T_FWHM_0 = (r2-r1) * 1e-12
T_e_0 = T_FWHM_0 / 2.35482


# Get FWHM of the output signal E
spline = UnivariateSpline(gp.timeax*1.0e12, np.abs(E[0]['E'][0])-1*np.max(np.abs(E[0]['E'][0]))/2, s=0)
r1, r2 = spline.roots() # find the roots
plt.annotate(s='', xy=(r1,np.max(np.abs(E[0]['E'][0]))/2), xytext=(r2,np.max(np.abs(E[0]['E'][0]))/2), arrowprops=dict(arrowstyle='<->'))
plt.text(r1+(r2-r1)/2.0, -0.04 + np.max(np.abs(E[0]['E'][0]))/2, '$T_{FWHM_1}$ = ' + str(np.round(r2-r1,2)) + ' ps', fontsize=12, horizontalalignment='center')
T_FWHM_1 = (r2-r1) * 1e-12
T_e_1 = T_FWHM_1 / 2.35482


# Get FWHM of the output signal E
spline = UnivariateSpline(gp.timeax*1.0e12, np.abs(Eq[0]['E'][0])-1*np.max(np.abs(Eq[0]['E'][0]))/2, s=0)
r1, r2 = spline.roots() # find the roots
plt.annotate(s='', xy=(r1,np.max(np.abs(Eq[0]['E'][0]))/2), xytext=(r2,np.max(np.abs(Eq[0]['E'][0]))/2), arrowprops=dict(arrowstyle='<->'))
plt.text(r1+(r2-r1)/2.0, -0.1 + np.max(np.abs(Eq[0]['E'][0]))/2, '$T_{FWHM_q}$ = ' + str(np.round(r2-r1,2)) + ' ps', fontsize=12, horizontalalignment='center')
T_FWHM_q = (r2-r1) * 1e-12
T_e_q = T_FWHM_q / 2.35482

plt.ylabel('$|E|$ a.u.'); plt.xlabel('Time $t$ [ps]'); legend = plt.legend(loc='upper right')

# Print the results
print 'Power: Gauss : T_FWHM_0 = ' , T_FWHM_0*1e12, ' ps'
print 'Power: Gauss : T_e_0 = ' , T_e_0*1e12, ' ps'
print 'Power: sqrt.Gauss : T_FWHM_1,0 = ' , T_FWHM_1*1e12, ' ps'
print 'Power: sqrt.Gauss : T_e_0,1 = ' , T_e_1*1e12, ' ps'
print 'Power: Gauss**2 : T_FWHM_0 = ' , T_FWHM_q*1e12, ' ps'
print 'Power: Gauss**2 : T_e_0 = ' , T_e_q*1e12, ' ps'

