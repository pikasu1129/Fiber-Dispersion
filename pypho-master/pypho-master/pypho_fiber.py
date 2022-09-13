#!/usr/bin/env python 
#  
#  Copyright 2014 Arne Striegler (arne.striegler@hm.edu)
#   
#   
#  
# Simulates a fiber
# uiukhuh
#
########################################################################

import numpy as np
import sys
from pypho_functions import *
from pypho_fiber_birefringence import pypho_fiber_birefringence
import time
import pyfftw
import copy

#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from itertools import product, combinations  

#import pylab as p
#import mpl_toolkits.mplot3d.axes3d as p3

########################################################################


class pypho_fiber(object):
    def __init__(self, glova = None, fibertype = None, D = None, S = None, gamma = None, alpha = None, l = None, useYPol = True, birefarray = None, phi_max=None):

        if glova == None:            
            print ("ERROR: You must define the global variables")
            sys.exit("PyPho stopped!")

        self.glova     = glova
        self.D         = None
        self.S         = None
        self.l         = None
        self.gamma     = None
        self.alpha     = None
        self.birefarray = None
        self.phi_max   = None
        self.useYPol   = True
        self.set(fibertype, D, S, gamma, alpha, l, useYPol, birefarray, phi_max)


########################################################################

    def __call__(self, E = None, fibertype = None, D = None, S = None, gamma = None, alpha = None, l = None, useYPol = True, birefarray = None, phi_max=None):

        self.set(fibertype, D, S, gamma, alpha, l, useYPol, birefarray)

        return self.transmit(E, fibertype, D, S, gamma, alpha, l, phi_max)

########################################################################
    
    def transmit(self, E = None, fibertype = None, D = None, S = None, gamma = None, alpha = None, l = None, phi_max=None):
        """Transmit the signal"""

        if E == None:
            print ("ERROR: You must define an optical signal")
            sys.exit("PyPho stopped!")

        if type(E) != list:
            E = [E]


        self.input_array_fx    = np.zeros(self.glova.sps*self.glova.nos, dtype=np.complex128)
        self.output_array_fx   = np.zeros(self.glova.sps*self.glova.nos, dtype=np.complex128)
        self.input_array_bx    = np.zeros(self.glova.sps*self.glova.nos, dtype=np.complex128)
        self.output_array_bx   = np.zeros(self.glova.sps*self.glova.nos, dtype=np.complex128)
        pyfftw.import_wisdom(self.glova.wisdom_pyfftw)
        self.fft_fwdx = pyfftw.FFTW(self.input_array_fx, self.output_array_fx, direction='FFTW_FORWARD', flags=[self.glova.fftw3_plan], threads=self.glova.fftw3_threads)
        self.fft_bwdx = pyfftw.FFTW(self.input_array_bx, self.output_array_bx, direction='FFTW_BACKWARD', flags=[self.glova.fftw3_plan], threads=self.glova.fftw3_threads)


        if self.useYPol:
            self.input_array_fy    = np.zeros(self.glova.sps*self.glova.nos, dtype=np.complex128)
            self.output_array_fy   = np.zeros(self.glova.sps*self.glova.nos, dtype=np.complex128)
            self.input_array_by    = np.zeros(self.glova.sps*self.glova.nos, dtype=np.complex128)
            self.output_array_by   = np.zeros(self.glova.sps*self.glova.nos, dtype=np.complex128)
            self.fft_fwdy = pyfftw.FFTW(self.input_array_fy, self.output_array_fy, direction='FFTW_FORWARD', flags=[self.glova.fftw3_plan], threads=self.glova.fftw3_threads)
            self.fft_bwdy = pyfftw.FFTW(self.input_array_by, self.output_array_by, direction='FFTW_BACKWARD', flags=[self.glova.fftw3_plan], threads=self.glova.fftw3_threads)


        self.set(fibertype, D, S, gamma, alpha, l, self.useYPol)
        
        self.E = E
        
        self.gamma_intern = self.gamma * 1e-3
        self.max_step = 100

        if self.useYPol:
            self.delta_beta = self.birefarray[0].delta_beta
            self.E[0]['E'] = self.coordrot(self.E[0]['E'], self.birefarray[0].angle)       
            self.birefindex = 1
            
        z = 0

        
        self.fibtrans(z)


        return self.E

     
 ########################################################################

    def fibtrans(self, zsep):
        """Calculate step"""
        
        lc = 0 
        lcd = 0 
       
        tic = time.time()

        self.E[zsep]['noise'] *= np.exp(-self.l * self.alpha)
        
        while lc < self.l:
            
            power_x = self.get_power_Ax(self.E[zsep]['E'])
            power_y = self.get_power_Ay(self.E[zsep]['E'])
            power = power_x + power_y
            # calculate the steplength

            next_leff = self.phi_max / np.max(power)         #(1-np.exp(-self.alpha*5))/self.alpha
            
            if next_leff*self.alpha < 1 :
                next_l = -np.log(1-next_leff*self.alpha)/self.alpha
            else:
                next_l = next_leff
                
            if next_l < next_leff:
                next_l = next_leff
            
            if next_l > self.max_step:
                next_l =  self.max_step
                next_leff = (1-np.exp(-self.alpha*next_l)) / self.alpha
            
            if lc + next_l > self.l:
                next_l =  self.l - lc
                next_leff = (1-np.exp(-self.alpha*next_l)) / self.alpha
                
            if self.useYPol:
                            
                doCoordRot = False
                if self.birefindex < len(self.birefarray) and lc + next_l >= self.birefarray[self.birefindex].z_point:
                    next_length = self.birefarray[self.birefindex].z_point - lc
                    self.birefindex += 1
                    next_l = next_length
                    next_leff = (1-np.exp(-self.alpha*next_l)) / self.alpha
                    doCoordRot = True

            # linear & nonlinear calculation
            #next_l = 100
            #next_leff = next_l            
            #self.E[zsep]['E'][0], self.E[zsep]['E'][1] = self.timestep(self.E[zsep]['E'][0], self.E[zsep]['E'][1], power_x, power_y, next_leff)
            
            self.E[zsep]['E'][0], self.E[zsep]['E'][1] = (self.E[zsep]['E'][0] * np.exp(-1j*next_leff*self.gamma_intern*(power_x + 0.6666666666666666 * power_y)), 
                self.E[zsep]['E'][1] * np.exp(-1j*next_leff*self.gamma_intern*(power_y + 0.6666666666666666 * power_x))) 

            self.E[zsep]['E'] = self.freqstep(self.E[zsep]['E'], next_l)
            

            if self.useYPol and doCoordRot:
                self.E[zsep]['E'] = self.coordrot(self.E[zsep]['E'], self.birefarray[self.birefindex-1].angle)
                self.delta_beta = self.birefarray[self.birefindex-1].delta_beta
                
                
            lc += next_l
            lcd += next_l
           # print (lc, next_l, next_leff)
            if lcd > 1000 :
                print (lc)
                lcd = 0
                
        if self.useYPol:
            total_Angle = 0
            for point in self.birefarray:
                total_Angle += point.angle
    
            self.E[zsep]['E'] = self.coordrot(self.E[zsep]['E'], -total_Angle)
               

        print (time.time() - tic )
        print ('OFF', lc)


########################################################################

    def freqstep(self, E, dz):
        """Calculate the linear step"""

        res = copy.deepcopy(E)

        self.input_array_fx[:] = E[0]
         
        pyfftw.FFTW.execute(self.fft_fwdx)
        
        self.input_array_bx[:] = self.output_array_fx * np.exp( self.beta_fac * dz) / len(E[0])
        
        if self.useYPol:
            self.input_array_bx[:] *= np.exp(-0.5 * 1j * fftshift (self.delta_beta * self.Domega) * dz)
        
        pyfftw.FFTW.execute(self.fft_bwdx)
        
        res[0] = self.output_array_bx * np.sqrt( np.exp(-self.alpha * dz) )

        if self.useYPol:
            self.input_array_fy[:] = E[1]
            
            
            pyfftw.FFTW.execute(self.fft_fwdy)
            
            self.input_array_by[:] = self.output_array_fy * np.exp( self.beta_fac * dz) / len(E[1])
            self.input_array_by[:] *= np.exp(0.5 * 1j * fftshift (self.delta_beta * self.Domega) * dz)
            
            pyfftw.FFTW.execute(self.fft_bwdy)
            
            res[1] = self.output_array_by * np.sqrt( np.exp(-self.alpha * dz) )

        return res
             

########################################################################

    def coordrot(self, E, angle):
        A = copy.deepcopy(E)
        A[0] = E[0]*np.cos(angle) - E[1]*np.sin(angle)
        A[1] = E[0]*np.sin(angle) + E[1]*np.cos(angle)
        return A

########################################################################

    def get_power_t(self, E):
        
        ox_real = np.real(E[0,:])
        ox_imag = np.imag(E[0,:])
        oy_real = np.real(E[1,:])
        oy_imag = np.imag(E[1,:])

        return ox_real**2 + ox_imag**2 + oy_real**2 + oy_imag**2

########################################################################
  
    def get_power_Ax(self, E):

        ox_real = np.real(E[0])
        ox_imag = np.imag(E[0])
        return  ox_real**2 + ox_imag**2

########################################################################

    def get_power_Ay(self, E):
        
        oy_real = np.real(E[1])
        oy_imag = np.imag(E[1])
        return oy_real**2 + oy_imag**2
               
########################################################################
        
    def get_parameters(self):
        """Get all fiber parameters set by the user"""
        biref_angle = []
        biref_zpoint = []
        biref_delta_beta = []
        for x in self.birefarray:
            biref_angle.append(x.angle)
            biref_zpoint.append(x.z_point)
            biref_delta_beta.append(x.delta_beta)

        params = {'type' : 'fiber', 'D':  self.D, 'S' : self.S,'l' : self.l, 'gamma' : self.gamma, 'alpha' : 4.343e3*self.alpha, 'biref_delta_beta' : biref_delta_beta, 'biref_zpoint' : biref_zpoint, 'biref_angle' : biref_angle, 'useYPol' : self.useYPol}
        return params

########################################################################
    def set(self, fibertype = None, D = None, S = None, gamma = None, alpha = None, l = None, useYPol = None, birefarray = None, phi_max=None):
        """Set fibre properties"""


        if fibertype ==  'SSMF':
            self.D         = 17.0                   # [ps/(nm km)]
            self.S         = 0.0                         # [ps/(nm2 km)]
            self.gamma     = 1.14                   # [1/(W m)]
            self.alpha     = db2neper(0.20)  / 1.0e3     # [1/m]


        if fibertype == 'DCF':
            self.D         = -100.0                    # [ps/(nm km)]
            self.S         = 0.0                        # [ps/(nm2 km)]
            self.gamma     = 1.7                        # [1/(W m)]
            self.alpha     = db2neper(0.20)  / 1.0e3   # [1/m]

        if D == None and fibertype == None and self.D == None:
            self.D             = 17.0
            print ("Warning: D and fibertype not specified, so I set D = ",self.D," ps/nm/km")
        elif D != None:
            self.D = D

        if S == None and fibertype == None and self.S == None:
            self.S             = 0.68
            print ("Warning: S and fibertype not specified, so I set S = ",self.S," ps/nm**2/km")            
        elif S != None:
            self.S = S

        if gamma == None and fibertype == None and self.gamma == None:
            self.gamma            = 1.14
            print ("Warning: gamma and fibertype not specified, so I set gamma = ",self.gamma," 1/W/km")            
        elif gamma != None:
            self.gamma = gamma

        if alpha == None and fibertype == None and self.alpha == None:
            print ("Warning: alpha and fibertype not specified, so I set alpha = 0.25 dB/km")
            self.alpha         = db2neper(0.25)  / 1.0e3
        elif alpha != None:
            self.alpha = db2neper(alpha) / 1.0e3
            
        if l == None and fibertype == None and self.l == None:
            self.l         = 80.0e3
            print ("Warning: L and fibertype not specified, so I set L = ",self.l," km")            
        elif l != None:
            self.l = l

        if useYPol == None and self.useYPol == None:
            self.useYPol = True
        elif useYPol != None:
            self.useYPol = useYPol
            
        if phi_max == None and self.phi_max == None:
            self.phi_max = 0.01
            print ("Warning: phi_max not specified, so I set phi_max = ", self.phi_max)
        elif phi_max != None:
            self.phi_max = phi_max            
        
        if birefarray == None and self.birefarray == None:
            nosiref = pypho_fiber_birefringence (0, 0, 0)
            self.birefarray = [nosiref]
        elif birefarray != None:
            self.birefarray = birefarray

        (self.beta_2, self.beta_3) =  DS2beta(self.D, self.S, self.glova.lambda0)
        self.Domega     = 2.0 * np.pi * (self.glova.freqax() - self.glova.f0) / 1.0e12
        #self.beta_fac   = -1j * fftshift  ( self.beta_2*0.5e24 * self.Domega**2.0 + self.beta_3*1.0e36 / 6.0  * self.Domega**3.0)
        Domega_tmp     = self.Domega**2.0 
        self.beta_fac       = self.beta_2*0.5e24 * Domega_tmp
        Domega_tmp     *= self.Domega
        self.beta_fac       += self.beta_3*1.0e36 / 6.0 * Domega_tmp
        self.beta_fac       = -1j * fftshift  (self.beta_fac)
        Domega_tmp     = 0
        
        y_max = self.glova.frange * 2.0 * np.pi  / 1.0e12
        x = np.arange(0,self.glova.sps*self.glova.nos)
        k =  2.0 * np.pi *self.glova.fres  / 1.0e12
        #self.y = x*k - ( (x*k-y_max/2.0)/np.abs(x*k-y_max/2.0) +1)*y_max/2.0
        
        

########################################################################



