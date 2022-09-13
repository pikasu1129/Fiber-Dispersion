#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho_arbmod.py
#
#  Copyright 2018 Arne Striegler (arne.striegler@hm.edu)
#
#
#
# Modulates in arbitray constellation points
#
#
########################################################################

import numpy as np
import sys
from pypho_functions import *

########################################################################

class pypho_arbmod(object):

    def __init__(self, glova = None):

        if glova == None:
            print ("ERROR: pypho_armod: You must define the global variables")
            sys.exit("PyPho stopped!")     
        
        self.glova         = glova


########################################################################

    def __call__(self,  glova = None, E = None, constpoints = None, symbols = None):

        if E == None:
            print ("ERROR: pypho_armod: You must define an optical signal E")
            sys.exit("PyPho stopped!")

        self.set(E, constpoints, symbols)
                
        for c in range(0, self.glova.nos):              

            E[0]['E'][0][c*self.glova.sps:(c+1)*self.glova.sps] *= self.constpoints[0][0][ int(symbols[0][c]) ]

            E[0]['E'][1][c*self.glova.sps:(c+1)*self.glova.sps] *= self.constpoints[1][0][ int(symbols[1][c]) ]
           
        return E

########################################################################

    def set(self, E = [], constpoints = None, symbols = None):


        if constpoints == None and self.constpoints == None:
            print ("ERROR: pypho_armod: No constellation points defined")
            sys.exit("PyPho stopped!")   
        elif constpoints != None:
            self.constpoints = constpoints
            self.M = int(np.floor( np.log2(len(constpoints[0][0])) ))

        if symbols == None and self.symbols == None:
            print ("ERROR: pypho_armod: No symbols sequence definded")
            sys.exit("PyPho stopped!")
        elif symbols != None:
            if len(symbols) == 1:
                symbols = [symbols, symbols]
                
            self.symbols = symbols
            

########################################################################
