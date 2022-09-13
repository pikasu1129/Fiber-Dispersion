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
from pypho_functions import *
import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

plt.close('all')

# Define network elements
gp       = pypho_setup(nos = 1*128, sps = 128, symbolrate = 10e9)
symsrc   = pypho_symbols(glova = gp, nos = gp.nos, pattern = 'random')
esigsrc  = pypho_signalsrc(glova = gp, pulseshape = 'gauss_rz' , fwhm = 0.33)
sig_1550 = pypho_lasmod(glova = gp, power = 0, Df = 0, teta = np.pi/4.0)
SSMF     = pypho_fiber(glova = gp, l = 10.0e3,  D = 17.0,   S = 0, alpha = 0.2, gamma = 1.14, phi_max = 0.1)
DCF      = pypho_fiber(glova = gp, l = 1.0e3,  D = -SSMF.D*SSMF.l*1.0e-3,   S = 0, alpha = 0.2e-12, gamma = 1.0e-9, phi_max = 10.0)
amp      = pypho_oamp(glova = gp, Pmean = 3.0, NF = 5)
osnr     = pypho_osnr(glova = gp)
modulator= pypho_arbmod(glova = gp) 

# Simulation

# Create symbolpattern
symbols_x = symsrc(p1=16)
symbols_y = symsrc()

# Create pulsetrain
onebits = symsrc(pattern = 'ones')
esig = esigsrc(bitsequence = onebits)
E_Tx = sig_1550(esig = esig)                                                      



# OOK
constpts_ook = [(  [0.001, 1.0]),
                (  [(0)], [(1)] ) ] 


# 8-PSK
constpts_8psk = [(  [ 1.0*np.exp(2.0j*np.pi*x/8.0) for x in range(0,8)] ),
                 (  [(0),(0),(0)], [(0),(0),(1)], [(0),(1),(0)], [(0),(1),(1)], [(1),(1),(1)], [(1),(1),(0)], [(1),(0),(1)], [(1),(0),(0)])] 

# 16-PSK
constpts_16psk = [(      [ 1.0*np.exp(2.0j*np.pi*x/16.0) for x in range(0,16)] ),
             ([(0),(0),(0),(0)], [(0),(0),(0),(1)], [(0),(0),(1),(0)], [(0),(0),(1),(1)], [(0),(1),(0),(0)], [(0),(1),(0),(1)], [(0),(1),(1),(0)], [(0),(1),(1),(1)],
              [(1),(1),(1),(1)], [(1),(1),(1),(0)], [(1),(1),(0),(1)], [(1),(1),(0),(0)], [(1),(0),(1),(1)], [(1),(0),(1),(0)], [(1),(0),(0),(1)], [(1),(0),(0),(0)]
             )]    

# 4-QAM
constpts_4qam = [(      [ 1.0*np.exp(2.0j*np.pi*x/4.0) for x in range(0,4)] ),
                 (      [(0),(0)], [(0),(1)], [(1),(1)], [(1),(0)]             )]    


# 16-QAM
alpha = np.arctan(np.sqrt(1.0)/3.0)
constpts_16qam = [(           [np.sqrt(3.0**2 + 1.0)]*8 + [np.sqrt(2.0)]*4 + [np.sqrt(2*3.0**2)]*4),
             (          [2.0*np.pi*x/4.0+alpha for x in range(0,4)] + [2.0*np.pi*x/4.0+np.pi-alpha for x in range(0,4)] + [2.0*np.pi*x/4.0+np.pi/4 for x in range(0,8)] ),

             ([(0),(0),(0),(0)], [(0),(0),(0),(1)], [(0),(0),(1),(0)], [(0),(0),(1),(1)], [(0),(1),(0),(0)], [(0),(1),(0),(1)], [(0),(1),(1),(0)], [(0),(1),(1),(1)],
              [(1),(1),(1),(1)], [(1),(1),(1),(0)], [(1),(1),(0),(1)], [(1),(1),(0),(0)], [(1),(0),(1),(1)], [(1),(0),(1),(0)], [(1),(0),(0),(1)], [(1),(0),(0),(0)]
             )]  # codes not optimized!  

constpts = constpts_16psk
print(constpts)
E = modulator( E = E_Tx, constpoints = [constpts, constpts], symbols = [symbols_x, symbols_y] )          # Modulate

P0 = 4    
E = amp(E = E, Pmean = P0)
E = osnr( E = E, OSNR = 58.0 )          # Set initial OSNR to 58 dB

plt.figure(1)  
plt.subplot(2, 1, 1); plt.grid(True); plt.title("Input signal", loc='left');plt.grid(True);
E_samp = E[0]['E'][0][int(gp.nos/2)::gp.nos]; plt.plot(np.real(E_samp), np.imag(E_samp), 'ro')   # Plot constallation diagramme

E_Tx = copy.deepcopy(E)
n_span = 0

E = amp(E = E, Pmean = P0)

for c in range(0, n_span):  # Transmission fiber
    print('Span: ', c)
    fibres = pypho_fiber_birefringence.create_pmd_fibre(SSMF.l, 1.0e3, 0.00)
    E = SSMF(E = E, birefarray = fibres)      
    E = amp(E = E, Pmean = P0)

for c in range(0, n_span):  # Dispersion compensation
    E = DCF(E = E)  

plt.figure(1)  
plt.subplot(2, 1, 2); plt.grid(True); plt.title("Output signal", loc='left');plt.grid(True);
E_samp = E[0]['E'][0][int(gp.nos/2)::gp.nos]; plt.plot(np.real(E_samp), np.imag(E_samp), 'go')   # Plot constallation diagramme



# Plot power and phase of the signal
plt.figure(3)
plt.subplot(2, 1, 1)
plt.title("Input signal", loc='left')
plt.plot(gp.timeax()*1.0e12, np.abs(E[0]['E'][0])**2, 'r', label='$E_x(0, t)$')
plt.plot(gp.timeax()*1.0e12, np.abs(E[0]['E'][1])**2, 'g', label='$E_y(0, t)$')
plt.ylabel('$10log |E_{x,y}|^2$'); plt.xlabel('Time [ps]'); plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(gp.timeax()*1.0e12, np.angle(E[0]['E'][0]), 'r')
plt.plot(gp.timeax()*1.0e12, np.angle(E[0]['E'][1]), 'g')
plt.ylabel('$ \phi_{x,y} $'); plt.xlabel('Time [ps]');plt.grid(True); plt.legend(); plt.show()