#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho.py
#
#  Copyright 2014 Arne Striegler (arne.striegler@hm.edu)
#
#
#

# Class pypho_setup setup
# Global parameter definition
#
########################################################################

import numpy as np
import sys
import os
import pyfftw
import uuid

from pypho_constants import pypho_constants
cv = pypho_constants();

########################################################################

class pypho_setup(object):

    def __init__(self, nos = None,sps = None, f0 = None, symbolrate = None, cloud = None ):

        self.nos         = None    # number of symbols
        self.sps         = None    # samples per symbol
        self.f0          = None    # center frequency in 1/s        
        self.symbolrate  = None    # reference symbolrate
        self.cloud       = None
        self.wisdom_pyfftw = None  # FFTW wisdom for pyfftw wrapper

        self.pswd = "testuser";
        self.key  = "testpwd";

        self.fftw3_threads = 8
        self.fftw3_plan = 'FFTW_MEASURE' #  FFTW_ESTIMATE,  FFTW_MEASURE,  FFTW_PATIENT,  FFTW_EXHAUSTIVE
        self.fftw3_path_pyfftw = "../config/wisdom_pyfftw_" + str(nos * sps); #" + hex(uuid.getnode()) + "
        self.fftw3_path_fftw3  = "../config/wisdom_fftw3.txt"; #+ hex(uuid.getnode()) + "_"

        # Check and create FFTW wisdom

        if  os.path.isfile(self.fftw3_path_pyfftw+".npy"):
            print('Found wisdom!')
            self.wisdom_pyfftw = np.load(self.fftw3_path_pyfftw+".npy")
            pyfftw.import_wisdom(self.wisdom_pyfftw)

        input_array  = np.zeros(nos * sps, dtype=np.complex128)
        output_array = np.zeros(nos * sps, dtype=np.complex128)
        print('Checking / creating wisdom to speed up FFTW...')
        fft_fwdx = pyfftw.FFTW(input_array, output_array, direction='FFTW_FORWARD', flags=[self.fftw3_plan], threads=self.fftw3_threads)
        print('Checking / creating wisdom to speed up FFTW...')
        fft_fwdy = pyfftw.FFTW(input_array, output_array, direction='FFTW_BACKWARD', flags=[self.fftw3_plan], threads=self.fftw3_threads)
        np.save(self.fftw3_path_pyfftw, pyfftw.export_wisdom())
        self.wisdom_pyfftw = pyfftw.export_wisdom()

        del(input_array, output_array, fft_fwdx, fft_fwdy)

        self.set(nos, sps, f0, symbolrate, cloud)

########################################################################



# frange: Frequency range
    @property
    def frange(self):
        """Get frequency range sps * symbolrate"""
        return self.sps * self.symbolrate

# sampletime: Sample time
    @property
    def sampletime(self):
        """Get sampletime symbolrate / sps"""
        return 1.0/self.symbolrate / self.sps

# timeax: Time axis
    def timeax(self):
        """Get timeax"""
        return np.linspace(self.sampletime, 1.0/self.symbolrate * self.nos, self.sps * self.nos) - self.sampletime

# freqax: Frequency axis
    def freqax(self):
        """Get frequency axis"""
        return np.arange(- self.frange/2.0, self.frange/2.0, self.fres) + self.f0

# fres: Frequency resolution
    @property
    def fres(self):
        """Get frequency resolution"""
        return self.symbolrate / (self.nos)


########################################################################
    def set(self, nos = None, sps = None, f0 = None, symbolrate = None, cloud = None ):
        """Set  properties"""          
            
        if nos == None and self.nos == None:            
            self.nos = 256
            print ("WARNING: Number of symbols not specified, so I am using ",self.nos ," nos")
        elif nos != None:
            self.nos = nos
            
        if self.nos < 0 or self.nos & (self.nos-1) != 0:        
            print ("Error: nos value must be 2**n")
            sys.exit("PyPho stopped!")            
            
        if sps == None and self.sps == None:
            self.sps = 256
            print ("WARNING: Samples per symbol not specified, so I am using ",self.sps," sps")
        elif sps != None:
            self.sps = sps       
            
        if f0 == None and self.f0 == None:        
            self.f0 = 193414489032258.06
            print ("WARNING: Center frequency f0 not specified, so I am using ", self.f0, " Hz")
        elif f0 != None:
            self.f0 = f0   
            
        self.lambda0     = cv.lightspeed/ self.f0    # center wavelength in m         
        
        if symbolrate == None and self.symbolrate == None:
            self.symbolrate = 10e9
            print ("WARNING: Symbolrate not specified, so I am using ",self.symbolrate," Symbols / sec")

        elif symbolrate != None:
            self.symbolrate = symbolrate  
            
        if cloud == None and self.cloud == None:
            print ("WARNING: Cloud not specified, so I am using False")
            self.cloud = False
        elif cloud != None:
            self.cloud = cloud  
