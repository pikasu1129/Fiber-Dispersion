#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho.py
#
#  Copyright 2018 Arne Striegler (arne.striegler@hm.edu)
#
#
#
# Plots an optical eye diagramme
#
#
########################################################################

import numpy as np
import sys
from pypho_functions import *
from pypho_signalsrc import pypho_signalsrc
import matplotlib.pyplot as plt

########################################################################

class pypho_eye(object):

    def __init__(self, glova = None, polarisation = None, style = None):

        if glova == None:
            print ("ERROR: pypho_optfi: You must define the global variables")
            sys.exit("PyPho stopped!")
        self.glova 	   = glova        
        self.polarisation = None
        self.figure 	   = None
        self.style 	   = None

        self.set( polarisation, style)

########################################################################

    def __call__(self,  E = None, polarisation = None, style = None):


        if len(E) != 2:
            print ("ERROR: You must define an optical signal")
            sys.exit("PyPho stopped!")

        self.set(polarisation, style)

        return self.out(E)


########################################################################

    def out(self, E = []):

        styles = self.style.split(',')
        if len(styles) == 1:
            styles.append(styles[0])
            
        if (self.polarisation == 'x,y') or (self.polarisation == 'y,x') :
            plotwhatlist = ['x', 'y']
        else:
            plotwhatlist = [self.polarisation]


        sc = 0
        for plotwhat in plotwhatlist:
            if plotwhat == 'x':
                signal = np.abs(E[0,:])**2
            elif plotwhat == 'y':
                signal = np.abs(E[1,:])**2
            elif (plotwhat == 'x+y') or (plotwhat == 'y+x'):
                signal = np.abs(E[0,:])**2 + np.abs(E[1,:])**2

            signal = np.reshape(signal, (self.glova.nos, self.glova.sps))*1e3

            for i in np.arange(1,self.glova.nos):
                plt.plot(self.glova.sampletime*1e12 	* np.transpose(np.arange(self.glova.sps)-self.glova.sps), np.transpose(signal[i,:]), styles[sc])
                plt.plot(self.glova.sampletime*1e12 	* np.transpose(np.arange(self.glova.sps)), np.transpose(signal[i,:]), styles[sc])
                plt.plot(self.glova.sampletime*1e12 	* np.transpose(np.arange(self.glova.sps)+self.glova.sps), np.transpose(signal[i,:]), styles[sc])

            plt.ylabel('Power [mW]'); plt.xlabel('Time [ps]')
            plt.grid(b=True, which='major', color='0.65', linestyle='--')
            sc +=1

########################################################################

    def set(self, polarisation = None, style = None):


        if polarisation == None and self.polarisation == None:
             self.polarisation = "x,y"
        elif polarisation != None:
            self.polarisation = polarisation

        if style == None and self.style == None:
             self.style = "r,g"
        elif style != None:
            self.style = style


########################################################################
