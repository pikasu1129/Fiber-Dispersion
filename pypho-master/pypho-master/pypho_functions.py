#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho.py
#
#  Copyright 2019 Arne Striegler (arne.striegler@hm.edu)
#
#
#
# Functions
#
#
########################################################################

import numpy as np
import scipy.fftpack
#import pyfftw

import json
import os
#import pycurl
#from StringIO import StringIO
from pypho_constants import pypho_constants
import zlib
import zipfile
cv = pypho_constants();
import copy 
import matplotlib.pyplot as plt

########################################################################
def dbm2w( power_dbm ):
   "Transfers power value from dBm to W"
   return  10**(power_dbm/10.0) / 1000.0

########################################################################
def w2dbm( power_w ):
   "Transfers power value from W to dBm"
   return 10.0*np.log10(power_w*1e3)

########################################################################
def polrot(sig, alpha):
    "Rotates the polarisation of the signal"
    sigout = np.cos(alpha)*sig[0,:] - np.sin(alpha)*sig[1,:]
    sigout = np.vstack((sigout, np.sin(alpha)*sig[0,:] + np.cos(alpha)*sig[1,:]))
    return sigout

########################################################################
def db2neper(alpha_dB):
    "Converts attenuation parameter from dB/km into 1/km"
    return alpha_dB/4.3429

########################################################################
def DS2beta(D, S, lambda_ref):
    "Converts D and S into beta_2 [s**2 / km] and beta_3 [s** 3 /km]"
    if ( D == 0.0) and (S == 0.0):
        return (0.0, 0.0)

    D *= 1.0e-6
    S *= 1.0e3
    fact = lambda_ref**2 / (2.0 * np.pi * cv.lightspeed)

    return (-fact*D, fact**2 * ( S + (2.0 * D / lambda_ref)))

########################################################################
def fft(E):
    "Calculate FFT"
    return scipy.fftpack.fft(E)

########################################################################
def ifft(E):
    "Calculate inverse FFT"
    return scipy.fftpack.ifft(E)

########################################################################
def fftshift(E):
    "Shift the zero-frequency component to the center of the spectrum"
    return scipy.fftpack.fftshift(E)

########################################################################
def supresslowpow(E):
    "Supress all array alements < x dBm"

    E[ np.where(np.log10(np.abs(E)**2)) < -10 ] = 1e-30 + 0*1j

    return E

########################################################################
def getpower_dBm(E):
    "Get mean power in dBm"
    return (10.0*np.log10(np.mean(np.abs(E[0,:])**2)), 10.0*np.log10(np.mean(np.abs(E[1,:])**2)) ) 

########################################################################
def getpower_W(E):
    "Get mean power in W"
    return (np.mean(np.abs(E[0,:])**2), np.mean(np.abs(E[1,:])**2))



########################################################################
def electrical_filter(glova, esig, B):
    'Get Gauß filterd electrical Signal'
    sig = []
    filfunc = np.exp(-((glova.freqax + (-glova.f0))/(B*2.0e9 /2.3582))**2/2.0)**2   
    sig = ifft((fft(esig)) * fftshift(filfunc))
        
    return sig

########################################################################
def get_decision_matrix(glova, E, constptsarray, symbols, ofil):
    'Create decision matrix'

    #Esiggi = copy.deepcopy(E)
    Esiggi = ofil(E = E )
    #I_x_I, I_x_Q, I_y_I, I_y_Q = ninetydeghybrid(glova, Esiggi, LO) # Get photo currents

    #Esiggi[0]['E'][0] = I_x_I +1.0j*I_x_Q        # Now as electrial signal
    #Esiggi[0]['E'][1] = I_y_I +1.0j*I_y_Q        # Now as electrial signal
    
    Esx = Esiggi[0]['E'][0][int(glova.sps/2)::glova.sps]
    Esy = Esiggi[0]['E'][1][int(glova.sps/2)::glova.sps]
    Nsx = symbols[0]
    Nsy = symbols[1]
    #print(Esx)
    del(symbols)
    N_raster    = 1000.0 # Anzahl der Rasterpunkte
 
    for c in [0, 1]:
        
        Dec         = np.zeros([int(N_raster), int(N_raster)])        
        Z_tmp       = np.zeros([int(N_raster), int(N_raster)])
        Z_max       = np.zeros([int(N_raster), int(N_raster)])
    
        if 0 == c: 
            Es = Esx
            Ns = Nsx
            constpts = constptsarray[0]
        else:
            Es = Esy
            Ns = Nsy
            constpts = constptsarray[1]
                    
        Es_maxima = np.sort(np.abs( np.array([np.abs(np.min(np.real(Es))), np.abs(np.max(np.real(Es))), np.abs(np.min(np.imag(Es))), np.abs(np.max(np.imag(Es)))]) ) )
        Es_re_min = -3*Es_maxima[-1]
        Es_re_max = +3*Es_maxima[-1] 
        Es_im_min = -3*Es_maxima[-1]
        Es_im_max = +3*Es_maxima[-1] 
        
        Es_re_ax = np.arange(Es_re_min, Es_re_max, (Es_re_max - Es_re_min) / N_raster )
        Es_re_ax = Es_re_ax[0:int(N_raster)]  # to be shure to have the correct number of cells
        Es_im_ax = np.arange(Es_im_min, Es_im_max, (Es_im_max - Es_im_min) / N_raster )    
        Es_im_ax = Es_im_ax[0:int(N_raster)]  # to be shure to have the correct number of cells
        xx, yy = np.meshgrid(Es_re_ax, Es_im_ax, sparse = True)    

        for const_num in range(0, len(constpts[0])):
            Z_tmp       = np.zeros([int(N_raster), int(N_raster)])
        
            for n_samp in range(0, len(Es)):
                
                if (const_num == Ns[n_samp]):
                    Z_tmp += np.exp( -( (np.real(Es[n_samp])-xx)**2 + (np.imag(Es[n_samp])-yy)**2 ) / 0.050 )
            

            Z_tmp /= np.max(Z_tmp)
            Z_dec_tmp = np.where(Z_max < Z_tmp)    
            Dec[Z_dec_tmp] = const_num    
            Z_max[Z_dec_tmp] = Z_tmp[Z_dec_tmp]
        
        if 0 == c: 
            Dec_x       = copy.deepcopy(Dec)
            Esx_re_ax   = copy.deepcopy(Es_re_ax)
            Esx_im_ax   = copy.deepcopy(Es_im_ax)
        else:
            Dec_y       = copy.deepcopy(Dec)
            Esy_re_ax   = copy.deepcopy(Es_re_ax)
            Esy_im_ax   = copy.deepcopy(Es_im_ax)
        
    return Dec_x, Esx_re_ax, Esx_im_ax, Dec_y, Esy_re_ax, Esy_im_ax

########################################################################
def ninetydeghybrid(gp, Esig, LO):
    '90 deg hybrid'

    for pol in [1, 0]:
        E_1 = 0.5 * ( Esig[0]['E'][pol] +      LO[0]['E'][pol] )
        E_2 = 0.5 * ( Esig[0]['E'][pol] -      LO[0]['E'][pol] )
        E_3 = 0.5 * ( Esig[0]['E'][pol] + 1.0j*LO[0]['E'][pol] )
        E_4 = 0.5 * ( Esig[0]['E'][pol] - 1.0j*LO[0]['E'][pol] )

        # Electrical filter
        
        E_1 = electrical_filter(glova = gp, esig = np.abs(E_1)**2, B = gp.symbolrate*2.0e-9)
        E_2 = electrical_filter(glova = gp, esig = np.abs(E_2)**2, B = gp.symbolrate*2.0e-9)
        E_3 = electrical_filter(glova = gp, esig = np.abs(E_3)**2, B = gp.symbolrate*2.0e-9)
        E_4 = electrical_filter(glova = gp, esig = np.abs(E_4)**2, B = gp.symbolrate*2.0e-9)

                
        I_I =  E_1 - E_2
        I_Q =  E_3 - E_4
       
        
        if pol == 0:
            I_x_I = I_I
            I_x_Q = I_Q
        else:
            I_y_I = I_I
            I_y_Q = I_Q
            
    return I_x_I, I_x_Q, I_y_I, I_y_Q

########################################################################
    
def create_optnoise (gp, P_sig, OSNR):
    'Create noise vector'
    E_N = np.zeros((2, gp.sps*gp.nos)) + 0.0j
    P_N = P_sig/10.0**(OSNR/10.0) * (gp.frange/12.5e9)    
    noisesamples = np.random.randn(4, gp.sps*gp.nos) * np.sqrt(P_N/4)
    E_N[0] = noisesamples[0] + 1.0j*noisesamples[1]
    E_N[1] = noisesamples[2] + 1.0j*noisesamples[3]

    return E_N

########################################################################
    
def calc_BER (gp, E, constpts, OSNR, Dec_x, Dec_y, Esx_re_ax, Esx_im_ax, Esy_re_ax, Esy_im_ax, M, ofil, symbols):
    'Calculate BER value of given signal'
    
    Esiggi = copy.deepcopy(E)
    del(E)
    Nsx = symbols[0]
    Nsy = symbols[1]
    del(symbols)
    Esx_re_ax_max = Esx_re_ax[1] - Esx_re_ax[0]
    Esx_im_ax_max = Esx_im_ax[1] - Esx_im_ax[0]
    Esy_re_ax_max = Esy_re_ax[1] - Esy_re_ax[0]
    Esy_im_ax_max = Esy_im_ax[1] - Esy_im_ax[0]    
    
    BER=[0,0]
    
    for t in range(0, M):
        
        if (int(t/float(M)*100.0)- int((t-1)/float(M)*100.0)) > 0:
            print('Progress: ', t/float(M)*100.0, '%')
        
        #create noise vectors
        E_tmp   = copy.deepcopy(Esiggi)
        E_N     = create_optnoise (gp, np.mean(abs(E_tmp[0]['E'][0]**2) + abs(E_tmp[0]['E'][1]**2)), OSNR) 
        
        
        plt.figure(1)
        # add noise
        E_tmp[0]['E'][0] +=  E_N[0]
        E_tmp[0]['E'][1] +=  E_N[1]
        
        E_tmp = ofil( E = E_tmp )
    
        # detect signal with noise
        #I_x_I, I_x_Q, I_y_I, I_y_Q = ninetydeghybrid(gp, E_tmp, LO)
        #E_tmp[0]['E'][0] = I_x_I + 1.0j*I_x_Q
        #E_tmp[0]['E'][1] = I_y_I + 1.0j*I_y_Q
        Esx = E_tmp[0]['E'][0][int(gp.sps/2)::gp.sps]
        Esy = E_tmp[0]['E'][1][int(gp.sps/2)::gp.sps]
      
        
        #plt.figure(1)
        #plt.subplot(2, 1, 1); plt.plot(np.real(Esx), np.imag(Esx), '.')
        #plt.subplot(2, 1, 2); plt.plot(np.real(Esy), np.imag(Esy), '.')
      
        # Calc BER            
        Ex_x_pos = np.floor( (np.real(Esx) - Esx_re_ax[0]) / Esx_re_ax_max )
        Ex_y_pos = np.floor( (np.imag(Esx) - Esx_im_ax[0]) / Esx_im_ax_max )
        Ey_x_pos = np.floor( (np.real(Esy) - Esy_re_ax[0]) / Esy_re_ax_max )
        Ey_y_pos = np.floor( (np.imag(Esy) - Esy_im_ax[0]) / Esy_im_ax_max )
    
    
        for n_samp in range( 0, len(Esx) ) :              

            # X-Pol            
            if int(Ex_x_pos[n_samp]) < 1000 and  int(Ex_x_pos[n_samp]) >= 0 and int(Ex_y_pos[n_samp]) < 1000 and  int(Ex_y_pos[n_samp] )>= 0:
                
                const_num_ist = int( Dec_x[int(Ex_y_pos[n_samp]), int(Ex_x_pos[n_samp])])
                BER[0] += np.cumsum( np.ceil(np.add(constpts[0][1][const_num_ist], constpts[0][1][int(Nsx[n_samp])] ) % 2.0) )[-1]  
                if t == 0:
                    if np.cumsum( np.ceil(np.add(constpts[0][1][int(const_num_ist)] ,constpts[0][1][int(Nsx[n_samp])] ) % 2.0) )[-1] > 0:
                        #print(0, constpts[0][4][int(const_num_ist)] ,constpts[0][4][int(Nsx[n_samp])], BER[1], const_num_ist, Nsx[n_samp], n_samp)

                        plt.subplot(2, 1, 1); plt.plot(np.real(Esx[n_samp]), np.imag(Esx[n_samp]), color=str(1-Nsx[n_samp]*1/len(constpts[0][1])), marker='x')
                    else:
                        plt.subplot(2, 1, 1); plt.plot(np.real(Esx[n_samp]), np.imag(Esx[n_samp]), color=str(1-Nsx[n_samp]*1/len(constpts[0][1])), marker='o', markersize=2)
            
            
            # Y-Pol
            if int(Ey_x_pos[n_samp]) < 1000 and  int(Ey_x_pos[n_samp]) >= 0 and int(Ey_y_pos[n_samp]) < 1000 and  int(Ey_y_pos[n_samp]) >= 0:
                
                const_num_ist = Dec_y[int(Ey_y_pos[n_samp]), int(Ey_x_pos[n_samp])]
                BER[1] += np.cumsum( np.ceil(np.add(constpts[1][1][int(const_num_ist)], constpts[1][1][int(Nsy[n_samp])] ) % 2.0) )[-1]
                if t == 0:
                    if np.cumsum( np.ceil(np.add(constpts[1][1][int(const_num_ist)] ,constpts[1][1][int(Nsy[n_samp])] ) % 2.0) )[-1] > 0:
                        #print(1, constpts[1][4][int(const_num_ist)] ,constpts[1][4][int(Nsy[n_samp])], BER[1], const_num_ist, Nsy[n_samp], n_samp)
                        plt.subplot(2, 1, 2); plt.plot(np.real(Esy[n_samp]), np.imag(Esy[n_samp]), color=str(1-Nsy[n_samp]*1/len(constpts[0][1])), marker='x')
                    else:
                        plt.subplot(2, 1, 2); plt.plot(np.real(Esy[n_samp]), np.imag(Esy[n_samp]), color=str(1-Nsy[n_samp]*1/len(constpts[0][1])), marker='o', markersize=2)
      
              
    BER[0] = np.log10(BER[0] / float(len(constpts[0][0]) * M * len(Esx)   ))
    BER[1] = np.log10(BER[1] / float(len(constpts[1][0]) * M * len(Esy)   ))
        
    return BER
