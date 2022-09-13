#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho.py
#  
#  Copyright 2014 Arne Striegler (arne.striegler@hm.edu)
#  
#   
#  
# Creates an electrical signal from bitpattern
#
#
########################################################################

import numpy as np
import sys

########################################################################

class pypho_signalsrc(object):
    def __init__(self, glova = None, pulseshape = [], fwhm = None):
            
        if glova == None:            
            print ("ERROR: You must define the global variables")
            sys.exit("PyPho stopped!")

        self.glova           = glova
        self.bitsequence     = []
        self.pulseshape      = None
        self.fwhm            = None
        
        self.set(pulseshape, fwhm)

########################################################################        
    def __call__(self, bitsequence = [], pulseshape = None, fwhm = None):        
    
        if len(bitsequence) == 0 :
            print ("Error: pypho_signalsrc: No bit pattern defined!")        
            sys.exit("PyPho stopped!")
        elif len(bitsequence) > 0:
            self.bitsequence = bitsequence
        
        self.set(pulseshape, fwhm)
    
     
        if self.pulseshape == "gauss_rz":
            self.tau = self.fwhm * 1/self.glova.symbolrate / 2.3582
            E = self.gauss_rz()        
                        
        elif self.pulseshape == "gauss_nrz":
            self.tau = self.fwhm * 1/self.glova.symbolrate / 2.3582
            E = self.gauss_nrz()    
                
        elif self.pulseshape == "sech_rz":
            self.tau = self.fwhm * 1/self.glova.symbolrate / 1.76
            E = self.sech_rz()    
   
        elif self.pulseshape == "rect":
            self.tau = 0
            E = self.rect()                
        else:
            print ("ERROR: No valid pulse shape")
            sys.exit("PyPho stopped!")        
        
        del self.bitsequence
        
        return E            
            

### Create Gauss RZ Pulse shape ########################################

    def gauss_rz(self):
        NoAb = 4
        refsig   = np.exp(-( ( self.glova.timeax()[int(self.glova.nos*self.glova.sps/2 - NoAb*self.glova.sps) : int(self.glova.nos*self.glova.sps/2 + NoAb*self.glova.sps)] - 1/self.glova.symbolrate*self.glova.nos/2 )/ self.tau)**2 / 2.0 )                
        gesig    = np.zeros(self.glova.sps*self.glova.nos + 3*2*NoAb*self.glova.sps)        
        
        bits     = np.append(self.bitsequence[-NoAb:], self.bitsequence)
        bits     = np.append(bits, self.bitsequence[0:NoAb])
        for ndx in range(0, len(bits)):
            if 1 == bits[ndx]:               
                gesig[int(self.glova.sps * (ndx)) : int( self.glova.sps * (ndx + 2*NoAb ) )] +=  refsig
        del bits
        del  refsig
        return gesig[int(2*NoAb*self.glova.sps - self.glova.sps/2): int(2*NoAb*self.glova.sps + self.glova.sps*self.glova.nos - self.glova.sps/2)]
    
    
### Create sech RZ Pulse shape ########################################

    def sech_rz(self):
        
        NoAb = 4
        refsig   = 1/np.cosh( (self.glova.timeax()[int(self.glova.nos*self.glova.sps/2 - NoAb*self.glova.sps) : int(self.glova.nos*self.glova.sps/2 + NoAb*self.glova.sps)] - 1/self.glova.symbolrate*self.glova.nos/2 )/ self.tau)       
        gesig    = np.zeros(self.glova.sps*self.glova.nos + 3*2*NoAb*self.glova.sps)
        
        bits     = np.append(self.bitsequence[-NoAb:], self.bitsequence)
        bits     = np.append(bits, self.bitsequence[0:NoAb])

        for ndx in range(0, len(bits)):
            if 1 == bits[ndx]:               
                gesig[self.glova.sps * (ndx) :  self.glova.sps * (ndx + 2*NoAb )] +=  refsig
        del bits
        del  refsig
        return gesig[int(2*NoAb*self.glova.sps - self.glova.sps/2) : int(2*NoAb*self.glova.sps + self.glova.sps*self.glova.nos - self.glova.sps/2)]     

### Create rectabgle pulse shape ########################################

    def rect(self):
        refrect  = np.concatenate((np.ones(self.glova.sps), np.zeros( self.glova.sps*(self.glova.nos -1 ) )), axis=0)
        
        gesig = 0*refrect
        
        for ndx in range(0, self.glova.nos ):
            if self.bitsequence[ndx] == 1:
                gesig += np.roll(refrect, self.glova.sps * ndx)
        return gesig    
        
### Create Gauss NRZ Pulse shape ########################################

    def gauss_nrz(self):
        refgaussinc        = np.exp(-( ( self.glova.timeax() - 1/self.glova.symbolrate*self.glova.nos )/ self.tau)**2 / 2.0 )
        refgaussinc        = np.roll(refgaussinc, int(self.glova.sps/2))
        refgaussinc[int(self.glova.sps/2) : int(self.glova.sps) ]  = 1
        
        refgaussdec        = np.exp(-(   self.glova.timeax() / self.tau)**2 / 2.0)
        refgaussdec        = np.roll(refgaussdec, int(self.glova.sps/2) )
        refgaussdec[0 : int(self.glova.sps/2) ]  = 1
        
        refgauss         = np.exp(-( ( self.glova.timeax() - 1/self.glova.symbolrate*self.glova.nos/2 )/ self.tau)**2 / 2.0 )
        refgauss         = np.roll(refgauss, int(self.glova.sps * self.glova.nos/2 + self.glova.sps/2) )
        
        refgausscomp    = np.roll(refgauss, int(- self.glova.sps/2) )
        refgausscomp[int(self.glova.sps/2+1):self.glova.sps]    = refgauss[1:int(self.glova.sps/2)]
        refgausscomp[int(self.glova.sps+1):]    = 0
        
        refzero            = np.zeros(self.glova.sps * self.glova.nos)
        refone            = np.concatenate([ np.ones(self.glova.sps) , np.zeros(self.glova.sps * (self.glova.nos-1)) ])

        
        gesig = 0*np.array(refgaussinc)
        bits_tmp = np.concatenate([ [0] , self.bitsequence])

        for x in range(1, self.glova.nos+1):

            if   bits_tmp[x-1]==0 and bits_tmp[x]==1:
                gesig += np.roll(refgaussinc, self.glova.sps * (x-1))
            elif bits_tmp[x-1]==1 and bits_tmp[x]==1:
                gesig += np.roll(refone, self.glova.sps * (x-1))
            elif bits_tmp[x-1]==1 and bits_tmp[x]==0:
                gesig += np.roll(refgaussdec, self.glova.sps * (x-1))
            elif bits_tmp[x-1]==0 and bits_tmp[x]==0:
                gesig += np.roll(refzero, self.glova.sps * (x-1))                

                
        return gesig    

########################################################################
    def set(self, pulseshape = None, fwhm = None):
        """Set  properties"""          
            
        if pulseshape == None and self.pulseshape == None:
            print ("WARNING: No pulseshape specified, so I am using 'gauss_rz'")
            self.pulseshape = "gauss_rz"
        elif pulseshape != None:
            self.pulseshape = pulseshape
            
            
        if fwhm == None and self.fwhm == None:
            print ("WARNING: No FWHM puls width specified, so I am using 0.3 of 1/symbolrate")
            fwhm = 0.333
        elif fwhm != None:
            self.fwhm = fwhm
            
