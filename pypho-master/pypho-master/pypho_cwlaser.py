#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho.py
#  
#  Copyright 2014 Arne Striegler (arne.striegler@hm.edu)
#  
#   
#  

# Generates a cw-light source
# 
#
########################################################################

import numpy as np
import sys
from pypho_functions import *

########################################################################

class pypho_cwlaser(object):
    def __init__(self, glova = None, power = None, teta = None, Df = None):
            
        if glova == None:            
            print ("ERROR: pypho_cwlaser: You must define the global variables")
            sys.exit("PyPho stopped!")
            
        self.glova = glova
        self.power     = None
        self.teta      = None
        self.Df        = None
        
        self.set(power, teta, Df)    
            

########################################################################

    def __call__(self, power = None, teta = None, Df = None):
            
        self.set(power, teta, Df)        
          
        if type(self.Df) != list:     
            self.Df = [self.Df]  

        E = [];
      
        for Dfi in self.Df:
            Etmp = {}
                        
            esig =  np.sqrt( dbm2w( self.power ) ) * np.exp(1j*2.0*np.pi*Dfi*1.0e9*self.glova.timeax())
            
            Etmp['E'] = np.array([np.zeros(self.glova.sps*self.glova.nos), np.zeros(self.glova.sps*self.glova.nos) ])  + 0j
            Etmp['E'][0][:] = esig*np.cos(self.teta) 
            Etmp['E'][1][:] = esig*np.sin(self.teta) 
            del esig
            Etmp['Df'] = Dfi
            Etmp['noise'] = np.zeros(self.glova.nos * self.glova.sps) + dbm2w( self.power ) / 10.0**(5.8) / 12.5e9 * self.glova.fres            
            E.append( Etmp )
            del Etmp

        return E     
                
            

########################################################################

    def set(self, power = None, teta = None, Df = None):
    
        if power == None and self.power == None:
            self.power = 0
            print ("WARNING: pypho_cwlaser: No mean power value specified, so I use P = ", self.power," dBm")
        elif power != None:
            self.power = power
                    
        if teta == None and self.teta == None:
            self.teta = 0
            print ("WARNING: pypho_cwlaser: No polarisation angle specified, so I use teta = ", self.teta)            
        elif teta != None:
            self.teta = teta        
                    
        if Df == None and self.Df == None:
            self.Df = 0
            print ("WARNING: pypho_cwlaser: No variation from center frequency specified, so I use Df = ", self.Df)            
        elif Df != None:
            self.Df = Df        
                            

########################################################################
