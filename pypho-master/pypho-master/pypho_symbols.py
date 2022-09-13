#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho.py
#  
#  Copyright 2014 Arne Striegler (arne.striegler@hm.edu)
#  
#  23.10.2019: Changed from pypho_bits to pypho_symbols
#  

# Create Symbol pattern
#
########################################################################
import numpy as np
import sys
import random
from pypho_functions import *

########################################################################

class pypho_symbols(object):
    def __init__(self, glova = None, nos = None, pattern = None, p1 = None, p2 = None):
        
        if glova == None:            
            print ("ERROR: You must define the global variables")
            sys.exit("PyPho stopped!")
            
        self.glova       = glova
        self.nos         = None
        self.pattern     = None
        self.p1          = None
        self.p2          = None
        self.set(nos, pattern, p1, p2)
        
########################################################################        
    def __call__(self, nos = None, pattern = None, p1 = None, p2 = None):
        print('------------');
        self.set(nos, pattern, p1, p2)
    
        if self.pattern == "singlepulse":
            symbols = np.zeros(self.nos, dtype=np.int)
            symbols[int(self.glova.nos/2)] = 1    
            
        elif self.pattern == "ones":
            symbols = np.ones(self.nos, dtype=np.int)
                        
        elif self.pattern == "random":    
            symbols         = self.rndm()       

        elif self.pattern == "debruijn":
            self.p2 = int( np.log10(self.glova.nos)/ np.log10(self.p1))
            print('Hint! I set Order n=p2=', self.p2) 
            if self.nos > 0:
                symbols = self.debruijn(self.p1, self.p2)[0:self.nos]
            else:
                symbols = self.debruijn(self.p1, self.p2)
             
        else:    
            print ("ERROR: No valid pattern specified!")
            sys.exit("PyPho stopped!")
            
    
        return symbols

### Create Random Bit Sequence #########################################

    def rndm(self):
        value_rndn = np.array( [[x]*int(np.ceil(self.glova.nos/self.p1)) for x in range(0, self.p1)] )
        value_rndn = value_rndn.reshape(-1)
        random.shuffle(value_rndn)
        value_rndn     = value_rndn[0:self.glova.nos]

        return value_rndn
    
### Create De Bruijn Sequence #########################################
  
    def debruijn(self, k, n):
        a = [0] * k * n
        sequence = []
        def db(t, p):
            if t > n:
                if n % p == 0:
                    for j in range(1, p + 1):
                        sequence.append(a[j])
            else:
                a[t] = a[t - p]
                db(t + 1, p)
                for j in range(a[t - p] + 1, k):
                    a[t] = j
                    db(t + 1, t)
        db(1, 1)
        return sequence    

########################################################################
    def set(self, nos = None, pattern = None, p1 = None, p2 = None):
        """Set  properties"""
        
        if nos == None and self.nos == None:
            print ("WARNING: Number of bit not specified, so I am using ", self.glova.nos)
            self.nos = self.glova.nos
        elif nos != None:
            self.nos = nos

        if pattern == None and self.pattern == None:
            self.pattern = 'random'
            print("WARNING: pattern not specified, so I am using ", self.pattern)
        elif pattern != None:
            self.pattern = pattern
            
        if p1 == None and self.p1 == None:
            self.p1 = 2
            print ("WARNING: p1 not specified, so I am using ", self.p1)
        elif p1 != None:
            self.p1 = p1            

        if p2 == None and self.p2 == None:
            self.p2 = int( np.log2(self.glova.nos) )
            print ("WARNING: p1 not specified, so I am using ", self.p2)
        elif p2 != None:
            self.p2 = p2
