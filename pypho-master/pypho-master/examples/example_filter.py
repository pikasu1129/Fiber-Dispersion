##!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Import functions and libraries
import sys
sys.path.append('../')
from pypho_setup import pypho_setup
from pypho_symbols import pypho_symbols
from pypho_signalsrc import pypho_signalsrc
from pypho_lasmod import pypho_lasmod
from pypho_oamp import pypho_oamp
from pypho_osnr import pypho_osnr
from pypho_functions import *
from pypho_optfi import *
import numpy as np
import copy
import matplotlib.pyplot as plt

from scipy.interpolate import UnivariateSpline

# Define network elements
gp       = pypho_setup(nos = 1*128, sps = 1*128, symbolrate = 10e9)
bitsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'random')
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.85)
sig_1550 = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = 0.25*np.pi)
osnr     = pypho_osnr(glova = gp)
opfilter = pypho_optfi(glova = gp, Df = 0, B = 50, filtype = 'cosrolloff', alpha = 0.8, loss = 0.0)
amp      = pypho_oamp(glova = gp, NF = 6)


# Simulation
bits    = bitsrc()
esig    = esigsrc(bitsequence = bits)
E       = sig_1550(esig = esig)                                                      

E = amp( E = E, Pmean = 0)          # Set mean power to 0 dBm
E = osnr( E = E, OSNR = 30.0 )      # Set initial OSNR to 30 dB
E_Tx = copy.deepcopy(E)

plt.figure(1)
plt.plot((gp.freqax()-gp.f0)*1.0e-9, 10.0*np.log10(E[0]['noise']*1e3), 'r', label='White gaussian noise')

E = copy.deepcopy(E_Tx)
E = opfilter(E = E, filtype = 'cosrolloff' )
plt.plot((gp.freqax()-gp.f0)*1.0e-9, 10.0*np.log10(E[0]['noise']*1e3), 'g', label='cosrolloff')

E = copy.deepcopy(E_Tx)
E = opfilter(E = E, filtype = 'gauss' )
plt.plot((gp.freqax()-gp.f0)*1.0e-9, 10.0*np.log10(E[0]['noise']*1e3), 'b', label='gauss')

E = copy.deepcopy(E_Tx)
E = opfilter(E = E, filtype = 'rect' )
plt.plot((gp.freqax()-gp.f0)*1.0e-9, 10.0*np.log10(E[0]['noise']*1e3), 'k', label='rect')

E = copy.deepcopy(E_Tx)
E = opfilter(E = E, filtype = 'gaussrolloff' )
plt.plot((gp.freqax()-gp.f0)*1.0e-9, 10.0*np.log10(E[0]['noise']*1e3), 'm', label='gaussrolloff')

E = copy.deepcopy(E_Tx)
E = opfilter(E = E, filtype = 'gauss_nth', alpha = 4 )
plt.plot((gp.freqax()-gp.f0)*1.0e-9, 10.0*np.log10(E[0]['noise']*1e3), 'c', label='gauss_nth')


plt.ylabel("Spectral power density [dBm/( " + str(gp.fres*1.0e-9) + " GHz)]"); plt.xlabel('Frequency deviation in GHz');
plt.grid(True)
plt.legend()
plt.show()