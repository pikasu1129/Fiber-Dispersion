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
from pypho_arbmod import pypho_arbmod
from pypho_oamp import pypho_oamp
from pypho_osnr import pypho_osnr
from pypho_optfi import pypho_optfi

from pypho_functions import *
import numpy as np
import copy
import matplotlib.pyplot as plt

from matplotlib import colors as mcolors
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

plt.close('all')

# Define network elements
gp          = pypho_setup(nos = 128, sps = 128, symbolrate = 10e9)
symsrc      = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'ones')
esigsrc     = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.25)
sig_1550    = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = np.pi/4.0)
SSMF        = pypho_fiber(glova = gp, l = 80.0e3,  D = 17.0,   S = 0, alpha = 0.3, gamma = 1.14, phi_max = 0.4)
DCF         = pypho_fiber(glova = gp, l = 1.0e3,  D = -SSMF.D*SSMF.l,   S = 0, alpha = 0.2e-12, gamma = 1.0e-12, phi_max = 10.0)
amp         = pypho_oamp(glova = gp, Pmean = 3.0, NF = 5)
osnr        = pypho_osnr(glova = gp)
modulator   = pypho_arbmod(glova = gp) 
filter_f0   = pypho_optfi(glova = gp, Df = 0, B = 50)

# Simulation

# Create symbol pattern
symbols_x = symsrc(pattern = 'random', p1 = 4)
symbols_y = symsrc(pattern = 'random', p1 = 4)

# Create pulsetrain
onebits = symsrc(pattern = 'ones')
esig    = esigsrc(bitsequence = onebits)
E_Tx    = sig_1550(esig = esig)                                                      

# Define constellation points: 4-QAM
constpts_x = [(      [ 1.0*np.exp(2.0j*np.pi*x/4.0+np.pi/4) for x in range(0, 4) ] ),
              (      [ (0),(0)], [(0),(1)], [(1),(1)], [(1),(0)] )]     

constpts_y = constpts_x         # Constellation 

E   = modulator( E = E_Tx, constpoints = [constpts_x, constpts_y], symbols = [symbols_x, symbols_y])          # Modulate

P0  = 0.0 
      
E   = amp(E = E, Pmean = 0)
E   = amp(E = E, Pmean = P0)
E   = osnr( E = E, OSNR = 10.0 )          # Set initial OSNR


for c in range(0, 0):
    print('Span: ', c)
    fibres = pypho_fiber_birefringence.create_pmd_fibre(SSMF.l, 1e3, 0.00)
    E = SSMF(E = E, birefarray = fibres)         
    E = DCF(E = E, l = 1.0)  
    E = amp(E = E, Pmean = P0)
    print('OSNR = ', osnr( E = E))


############################
# Calculate decision matrix
############################

E = amp(E = E, Pmean = 0)

Ein = copy.deepcopy(E)
plt.figure(1)
# Get decision matrix
Dec_x, Esx_re_ax, Esx_im_ax, Dec_y, Esy_re_ax, Esy_im_ax = get_decision_matrix(gp, Ein, [constpts_x, constpts_y], [symbols_x, symbols_y], filter_f0)

# Plot decision matrix
plt.figure(1)
plt.subplot(2, 1, 1); h = plt.contourf(Esx_re_ax, Esx_im_ax, Dec_x, 32, cmap=plt.cm.jet )
plt.plot(np.real(E[0]['E'][0][int(gp.sps/2)::gp.sps]), np.imag(E[0]['E'][0][int(gp.sps/2)::gp.sps]), '*')
plt.subplot(2, 1, 2); h = plt.contourf(Esy_re_ax, Esy_im_ax, Dec_y, 32, cmap=plt.cm.jet )
plt.plot(np.real(E[0]['E'][1][int(gp.sps/2)::gp.sps]), np.imag(E[0]['E'][1][int(gp.sps/2)::gp.sps]), '.')
############################
# Calculate BER
############################
BER = calc_BER (gp, E, [constpts_x, constpts_y], osnr( E = E), Dec_x, Dec_y, Esx_re_ax, Esx_im_ax, Esy_re_ax, Esy_im_ax, 10, filter_f0, [symbols_x, symbols_y])

print(BER)
plt.show()