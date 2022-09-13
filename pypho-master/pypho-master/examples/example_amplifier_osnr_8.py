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
from pypho_oamp import pypho_oamp
from pypho_osnr import pypho_osnr

import numpy as np
import matplotlib.pyplot as plt

# Define network elements
gp       = pypho_setup(nos = 128, sps = 128, symbolrate = 10e9)
bitsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'random')
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.85)
sig_1550 = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = 0.25*np.pi)
SSMF     = pypho_fiber(glova = gp, l = 8.0e3,  D = 17,   S = 0, alpha = 2.0, gamma = 1.0e-12, phi_max = 10.0) # High loss fiber, to speed up the simulation
amp      = pypho_oamp(glova = gp, NF = 6)
osnr     = pypho_osnr(glova = gp)

# Simulation
bits    = bitsrc()
esig    = esigsrc(bitsequence = bits)
E       = sig_1550(esig = esig)                                                      
E = osnr( E = E, OSNR = 58.0 )          # Set initial OSNR to 58 dB

P = 3.0                                 # Signal power

plt.figure(1); plt.plot(0, osnr(E = E ), 'r*')

for t in range(10):                     # Loop through 10 spans
    E = SSMF(E = E) 
    E = amp( E = E, Pmean = P)          # Set mean power to 3 dBm
    plt.plot(t+1, osnr( E = E ), 'r*')

plt.title("OSNR as function of number of spans", loc='left');plt.ylabel('OSNR [dB]');plt.xlabel('Span number'); plt.grid(True); plt.show()