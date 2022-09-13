# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 13:44:16 2018

@author: lonstud7
"""

import numpy as np
import sys
from pypho_functions import *
from pypho_constants import pypho_constants
cv = pypho_constants();


########################################################################

class pypho_oamp(object):
    def __init__(self, glova = None, G = None, Pmean = None, Pmax = None, NF = None):

        if glova == None:
            print ("ERROR: pypho_amp: You must define the global variables")
            sys.exit("PyPho stopped!")
        else:
            self.glova = glova

        self.G      = None
        self.Pmean  = None
        self.Pmax   = None
        self.NF     = None
        self.mode = "none" 
        
        self.set(G, Pmean, Pmax, NF)



########################################################################
    def __call__(self, E = None, G = None, Pmean = None, Pmax = None, NF = None):

        self.set(G, Pmean, Pmax, NF)

        return self.out(E = E)

########################################################################
    def out(self, E = None):

        if E == None:
            print ("ERROR: pypho_amp: You must define an optical signal")
            sys.exit("PyPho stopped!")

        if type(E) != list:
            E = [E]

        z = 0

        for Ei in E:
            if self.mode == "Pmean":
                (P1, P2) = getpower_W(E[z]['E'])
                self.G = 10.0 *np.log10( 1.0e-3*10.0**(self.Pmean/10.0) /  (P1+P2) )
                
            elif self.mode == "Pmax":
                P1 = np.max(np.abs(E[z]['E'][0])**2)
                P2 = np.max(np.abs(E[z]['E'][1])**2)
                self.G =  10.0 * np.log10( 1.0e-3*10.0**(self.Pmax / 10.0 ) / np.max( np.array([P1, P2]) ) ) 
            print ('G = ', self.G)
            E[z]['E']       = E[z]['E']     * np.sqrt( 10.0**(self.G/10.0) )
            E[z]['noise']   = E[z]['noise'] * 10.0**(self.G/10.0)

            if self.G > 0:
                E[z]['noise']   += ( 10.0**(self.G/10.0) - 1 ) * 10.0**(self.NF/10.0) * cv.pwquant * self.glova.f0 * self.glova.fres                
            else :
                print ("WARNING: pypho_oamp: Gain G < 0dB! Signal is only attenuated.")
       
            z += 1
        return E


########################################################################
    def set(self, G = None, Pmean = None, Pmax = None, NF = None):
        """Set amplifer parameter"""

        if NF == None and self.NF == None:
            self.NF = 4
            print ("WARNING: pypho_oamp: Noise figure NF not defined! I set NF =", self.NF, " dB")
        elif NF != None:
            self.NF = NF


        if G == None and Pmean == None and Pmax == None and self.G == None and self.Pmean == None and self.Pmax == None:
            self.G = 0
            print ("WARNING: pypho_oamp: No gain or power value defined! I set G =", self.G, " dB")
            self.mode = "G"

        if G != None :
            self.G = G
            self.mode = "G"
        elif Pmean != None:
            self.Pmean = Pmean
            self.mode = "Pmean"
        elif Pmax != None:
            self.Pmax = Pmax
            self.mode = "Pmax"

        if self.mode == "none":
            self.G = 0
            print ("WARNING!: pypho_oamp: No gain or power value defined! I set G =", self.G, " dB")
            self.mode = "G"
