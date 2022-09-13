#!/usr/bin/env python
#
#  Copyright 2018 Arne Striegler (arne.striegler@hm.edu)
#
#
#
# Simulates a fiber ...fast!
#
#
########################################################################
import numpy as np
import sys
from pypho_functions import *
from pypho_fiber_birefringence import pypho_fiber_birefringence
import time

import copy
from speedfiber import *
#from speedtest import *

########################################################################

class pypho_fiber(object):
    def __init__(self, glova = None, fibertype = None, D = None, S = None, gamma = None, alpha = None, l = None, birefarray = None, phi_max = None):

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
        self.nos        = None
        self.phi_max     = None

        self.set(fibertype, D, S, gamma, alpha, l, birefarray, phi_max)


########################################################################

    def __call__(self, E = None, fibertype = None, D = None, S = None, gamma = None, alpha = None, l = None, birefarray = None, phi_max = None):

        self.set(fibertype, D, S, gamma, alpha, l, birefarray, phi_max)

        if self.glova.cloud :
            E.append(self.get_parameters())
            return E
        else :
             
            if E == None:
                print ("ERROR: You must define an optical signal")
                sys.exit("PyPho stopped!")
    
            if type(E) != list:
                E = [E]
    
            self.set(fibertype, D, S, gamma, alpha, l, phi_max)
    
            #self.E = E
    
            self.gamma_intern = self.gamma * 1e-3
            self.max_step = 200
    
            z = 0
            

            tic0 = time.time()   
            n = self.glova.sps*self.glova.nos
            #Ef_out = np.zeros((1,n)) + 1j*np.ones((1,n))
            Ex_out = np.zeros(n) + 1j*np.ones(n)
            Ey_out = np.zeros(n) + 1j*np.ones(n)
            birefarraydoubles = np.zeros((len(self.birefarray), 3))
            for i in range (0, len(self.birefarray)):
                birefarraydoubles[i,0] = self.birefarray[i].angle
                birefarraydoubles[i,1] = self.birefarray[i].z_point
                birefarraydoubles[i,2] = self.birefarray[i].delta_beta
                
            cyfiber(self.glova.sps*self.glova.nos, self.l, np.asarray(E[z]['E'][0]), np.asarray(E[z]['E'][1]),
                    self.alpha, self.gamma_intern, self.phi_max, birefarraydoubles, len(self.birefarray), self.max_step, self.beta_2, self.beta_3,
                    self.glova.frange, self.glova.fres, self.glova.fftw3_threads, self.glova.fftw3_path_fftw3, Ex_out, Ey_out)
            #print(self.glova.fftw3_path_fftw3)
            E[z]['E'][0] = Ex_out
            E[z]['E'][1] = Ey_out
            E[z]['noise'] *= np.exp(-self.l * self.alpha)
            print ('Fertig: ', time.time() - tic0 )        
            del (Ex_out, Ey_out)
            return E

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

        params = {'type' : 'fiber', 'D':  self.D, 'S' : self.S,'l' : self.l, 'gamma' : self.gamma, 'alpha' : 4.343e3*self.alpha, 'biref_delta_beta' : biref_delta_beta, 'biref_zpoint' : biref_zpoint, 'biref_angle' : biref_angle, 'useYPol' : True}
        return params

########################################################################
    def set(self, fibertype = None, D = None, S = None, gamma = None, alpha = None, l = None, birefarray = None, phi_max = None):
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
            print ("Warning: D and fibertype not specified, so I set D = 16.8 ps/nm/km")
            self.D             = 17.0
        elif D != None:
            self.D = D

        if S == None and fibertype == None and self.S == None:
            print ("Warning: S and fibertype not specified, so I set S = 0.68 ps/nm**2/km")
            self.S             = 0.68
        elif S != None:
            self.S = S

        if gamma == None and fibertype == None and self.gamma == None:
            print ("Warning: gamma and fibertype not specified, so I set gamma = 1.14 1/W/km")
            self.gamma            = 1.14
        elif gamma != None:
            self.gamma = gamma

        if alpha == None and fibertype == None and self.alpha == None:
            print ("Warning: alpha and fibertype not specified, so I set alpha = 0.25 dB/km")
            self.alpha         = db2neper(0.25)  / 1.0e3
        elif alpha != None:
            self.alpha = db2neper(alpha) / 1.0e3

        if l == None and fibertype == None and self.l == None:
            print ("Warning: L and fibertype not specified, so I set L = 80 km")
            self.l         = 80.0e3
        elif l != None:
            self.l = l
            
        if phi_max == None and self.phi_max == None:
            print ("Warning: phi_max not specified, so I set phi_max = 1e-3")
            self.phi_max         = 1.0e-3
        elif phi_max != None:
            self.phi_max = phi_max

        if birefarray == None and self.birefarray == None:
            nosiref = pypho_fiber_birefringence (0.0, 0.0, 0.0)
            self.birefarray = [nosiref]
        elif birefarray != None:
            self.birefarray = birefarray

        (self.beta_2, self.beta_3) =  DS2beta(self.D, self.S, self.glova.lambda0)



########################################################################


