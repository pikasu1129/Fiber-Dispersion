#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho.py
#
#  Copyright 2018 Arne Striegler (arne.striegler@hm.edu)
#
#
#
#
#
# Simulations direct modulated laser source
# Without Chirp
#
########################################################################

import numpy as np
import sys
from pypho_functions import *
from pypho_constants import pypho_constants
cv = pypho_constants();

########################################################################

class pypho_lasmod(object):
    def __init__(self, glova = None, power = None, Df = None, teta = None):

        if glova == None:
            print ("ERROR: You must define the global variables")
            sys.exit("PyPho stopped!")

        self.glova     = glova
        self.power     = None
        self.esig     = []
        self.Df         = None
        self.teta     = None

        self.set(power, Df, teta)

########################################################################
    def __call__(self, esig = [], power = None, Df = None, teta = None):

        if len(esig) == 0 :
            print ("Error: pypho_lasmod: No electrical signal defined!")        
            sys.exit("PyPho stopped!")
        
        self.set(power, Df, teta)
       
        if type(self.Df) != list:     
            self.Df = [self.Df]  

        E = [];
        esig =  esig + 0.0j
        for Dfi in self.Df:
            Etmp = dict()
                            
            esig *= np.sqrt( dbm2w( self.power ) / np.mean(esig) ) * np.exp(1.0j*2.0*np.pi*Dfi*1.0e9*self.glova.timeax())
            
            Etmp['E'] = np.zeros((2, self.glova.nos * self.glova.sps)) + 0.0j
            Etmp['E'][0][:] = esig * np.cos(self.teta) 
            Etmp['E'][1][:] = esig * np.sin(self.teta) 
            del esig
            Etmp['Df'] = Dfi
            Etmp['noise'] = np.zeros(self.glova.nos * self.glova.sps) + dbm2w( self.power ) / 10.0**(5.8) / 12.5e9 * self.glova.fres            
            E.append( Etmp )
            del Etmp

        return E


########################################################################
    def set(self, power = None, Df = None, teta = None):
        """Set  properties"""
      
        if power == None and self.power == None:
            self.power = 0
            print ("WARNING: pypho_lasmod: No mean power value specified, so I use P = ",self.power ," dBm")
        elif power != None:
            self.power = power


        if teta == None and self.teta == None:
            self.teta = 0
            print ("WARNING: pypho_lasmod: No polarisation angle specified, so I use teta = ", self.teta)

        elif teta != None:
            self.teta = teta


        if Df == None and self.Df == None:
            print ("WARNING: pypho_cwlaser: No variation from center frequency specified, so I use Df = 0")
            self.Df = np.array([0])
        elif Df != None:
            self.Df = Df
