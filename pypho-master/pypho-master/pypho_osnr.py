# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 13:44:16 2018

@author: lonstud7
"""

import numpy as np
import sys
from pypho_functions import *
from pypho_constants import pypho_constants


########################################################################

class pypho_osnr(object):
    def __init__(self, glova = None, OSNR = None, Df = None):

        if glova == None:
            print ("ERROR: pypho_osnr: You must define the global variables")
            sys.exit("PyPho stopped!")
        else:
            self.glova = glova

        self.OSNR   = None
        self.Df      = None
        self.mode   = "get"
        self.glova  = glova
        
        self.set(OSNR, Df)

########################################################################
    def __call__(self, E = None, OSNR = None, Df = None):
        
        self.set(OSNR, Df)
        
        if E == None:
            print ("ERROR: pypho_osnr: You must define an optical signal")
            sys.exit("PyPho stopped!")

        if type(E) != list:
            E = [E]

        idx_f = (np.abs(self.glova.freqax() - self.Df)).argmin()
    
        z = 0
        for Ei in E:
            
            (P1, P2) = getpower_W(E[z]['E']) 

            if self.mode == "get":
                
                E = 10.0 * np.log10( (P1+P2) /  (E[z]['noise'][idx_f] * 12.5e9 / self.glova.fres ) )
                
            elif self.mode == "set":
                                
                k_noise = (P1+P2) / 10**(self.OSNR/10) / 12.5e9 * self.glova.fres
                
                E[z]['noise'] *=  k_noise / E[z]['noise'][idx_f]
                
            z += 1   
            
        return E


########################################################################
    def set(self, OSNR = None, Df = None):
        """Set or get OSNR value"""

        self.OSNR = OSNR
                
        if Df == None and self.Df == None:
            self.Df = 0
            print ("WARNING: pypho_osnr: Frequency deviation Df not defined! I set Df =", self.Df, " Hz")
        elif Df != None:            
            self.Df = Df
            
        if self.Df < -self.glova.frange/2.0 :
            self.Df = -self.glova.frange/2.0
            print ("WARNING: pypho_osnr: Frequency deviation Df smaller than min ! I set Df =", self.Df, " Hz")
        elif self.Df > self.glova.frange/2.0 :
            self.Df = self.glova.frange/2.0
            print ("WARNING: pypho_osnr: Frequency deviation Df  higher than max ! I set Df =", self.Df, " Hz")
        
        if OSNR == None:
            self.mode = "get"
        else:        
            self.mode = "set"