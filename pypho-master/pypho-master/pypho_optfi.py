#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pypho.py
#  
#  Copyright 2014 Arne Striegler (arne.striegler@hm.edu)
#  
#   
#  
# Simulations direct modulated laser source
# Without Chirp
#
########################################################################
import matplotlib.pyplot as plt
import numpy as np
import sys
from pypho_functions import *
import pyfftw

########################################################################

class pypho_optfi(object):
    def __init__(self, glova = None, Df = None, B = None, filtype = None, alpha = None, loss = None):

        if glova == None:            
            print ("ERROR: pypho_optfi: You must define the global variables")
            sys.exit("PyPho stopped!")
            
        self.glova      = glova
        self.Df         = None
        self.B          = None
        self.filtype    = None
        self.alpha      = None
        self.loss       = None

        self.set(Df, B, filtype, alpha, loss)                

########################################################################

    def __call__(self,  E = None, Df = None, B = None, filtype = None, alpha = None, loss = None):

        n_samples = self.glova.sps * self.glova.nos

        if E == None:
            print ("ERROR: pypho_optfi: You must define an optical signal")
            sys.exit("PyPho stopped!")
        
        
        self.set(Df, B, filtype, alpha, loss)    

        
        if self.filtype == "gaussrolloff":        
            z=0
            filfunc=np.zeros(self.glova.nos*self.glova.sps)
            res = (self.glova.symbolrate/self.glova.nos)
            offset = -self.Df*1e9-(self.glova.sps*self.glova.nos/2)*(self.glova.symbolrate/self.glova.nos)


            while z<(self.glova.nos*self.glova.sps):
                if(abs(z*res+offset)<=((1-self.alpha)*self.B/2*1e9)):
                    filfunc[z]=1
                    
                elif(0<z*res+offset-(1-self.alpha)*self.B/2*1e9):
                    filfunc[z]=np.exp(-(((z*res+offset)-(1-self.alpha)*self.B/2*1e9)*((np.sqrt(np.abs(np.log(0.5)))/(self.alpha*self.B/2*1e9))))**2)
            
                elif(0>z*res+offset+(1-self.alpha)*self.B/2*1e9):
                    filfunc[z]=np.exp(-(((z*res+offset)+(1-self.alpha)*self.B/2*1e9)*((np.sqrt(np.abs(np.log(0.5)))/(self.alpha*self.B/2*1e9))))**2)
                
                z=z+1        
        
        
            z = 0
            
        
        
        
        if self.filtype == "gauss":
            
            filfunc = np.exp(-( ( self.glova.freqax() + ( - self.glova.f0 - self.Df*1e9 ) )/ (self.B*2.0e9 / 2.3582))**2 / 2.0 )**2
            
            
            
        if self.filtype == "rect":
                        
            filfunc = np.zeros(self.glova.nos*self.glova.sps)+1e-12
            
            filfunc[ np.where(np.abs(self.glova.freqax()-self.glova.f0 -self.Df*1e9) <= self.B*0.5e9) ] = 1
            


        if self.filtype == "gauss_nth":
    
            filfunc = ( np.exp( -2.0*np.log(2.0) * ( (  - self.glova.freqax() + self.glova.f0 + self.Df*1e9 )  /  (self.B*0.5e9 / (2.0**(self.alpha-1) )**(1.0/float(self.alpha)) ) )**(2*self.alpha) ) )
 


        if self.filtype == "cosrolloff":
            
            T1      = 1.0 / (self.B * 1.0e9)

            offset  = -self.Df*1.0e9-(self.glova.sps*self.glova.nos/2.0)*(self.glova.symbolrate/self.glova.nos)

            z = 0
            
            filfunc = np.zeros(self.glova.nos*self.glova.sps)



            while z<(self.glova.nos*self.glova.sps):
                if(abs(z*self.glova.fres+offset)<=((1.0-self.alpha)/(2*T1))):
                    filfunc[z]=1.0
                elif (((1-self.alpha)/(2.0*T1)<abs(z*self.glova.fres+offset)) and (abs(z*self.glova.fres+offset))<=(1+self.alpha)/(2*T1)):
                    filfunc[z]=(np.cos((np.pi*T1)/(2.0*self.alpha)*(abs(z*self.glova.fres+offset)-((1-self.alpha)/(2*T1))))**2)
                else:
                    filfunc[z]= 0.0
    
                z += 1
            
        z = 0
        
        filfunc *= dbm2w(-1*self.loss)*1000.0
        #print("--->"+ self.glova.wisdom_pyfftw)
        pyfftw.import_wisdom(self.glova.wisdom_pyfftw)
        input_array_f  = np.zeros(n_samples, dtype=np.complex128)
        output_array_f = np.zeros(n_samples, dtype=np.complex128)
        fft_f = pyfftw.FFTW(input_array_f, output_array_f, direction='FFTW_FORWARD',  flags=[self.glova.fftw3_plan], threads=self.glova.fftw3_threads)
        input_array_b  = np.zeros(n_samples, dtype=np.complex128)
        output_array_b = np.zeros(n_samples, dtype=np.complex128)
        fft_b = pyfftw.FFTW(input_array_b, output_array_b, direction='FFTW_BACKWARD', flags=[self.glova.fftw3_plan], threads=self.glova.fftw3_threads)

        for Ei in E:
            input_array_f[:] = E[z]['E'][0]
            pyfftw.FFTW.execute(fft_f)
            input_array_b[:] = output_array_f * fftshift(np.sqrt( filfunc) )
            pyfftw.FFTW.execute(fft_b)
            E[z]['E'][0] = output_array_b[:] / n_samples


            input_array_f[:] = E[z]['E'][1]
            pyfftw.FFTW.execute(fft_f)
            input_array_b[:]  = output_array_f * fftshift(np.sqrt( filfunc) )
            pyfftw.FFTW.execute(fft_b)
            E[z]['E'][1] = output_array_b[:] / n_samples

           # E[z]['E'][0]    = ifft( ( fft(E[z]['E'][0])  )  *  fftshift(np.sqrt( filfunc) ))
           # E[z]['E'][1]    = ifft( ( fft(E[z]['E'][1])  )  *  fftshift(np.sqrt( filfunc) ))
            E[z]['noise']   *= filfunc
        
            z += 1  
          # TODO : Freemem FFTW
          # TODO : for Ei in E ändern in z 0 to length.E
          
        #fig=plt.figure('filter')
        #ax = fig.add_subplot(111)
        #ax.plot((self.glova.freqax()-self.glova.f0)*1e-9, filfunc)
        #ax.plot((self.glova.freqax() - self.glova.f0) * 1e-9, np.abs(input_array_f[:]))

        #ax.set_ylim(-0.2, 1.2)
        
        return E
            

########################################################################

    def set(self, Df = None, B = None, filtype = None, alpha = None, loss = None):
        
        
        if filtype == None and self.filtype != None:
            filtype = self.filtype
        
        if filtype == None and self.filtype == None:
            filtype = "gauss"
            print ("WARNING: pypho_optfi: No filter type specified! So I'm using gaussian fiter!")
   
        self.filtype = filtype            
        
        if filtype == "cosrolloff" or filtype == "gaussrolloff" or filtype == "gauss_nth":
            
            if alpha == None and self.alpha != None:
                alpha = self.alpha
            if alpha == None and self.alpha == None:
                print ("WARNING: pypho_optfi: Alpha is not specified! So I'm settin it to 0.25!")
                alpha = 0.25
                
            self.alpha = alpha
                
                   
        if Df == None and self.Df != None:            
            Df = self.Df
                    
        if Df == None and self.Df == None:            
            print ("WARNING: pypho_optfi: No deviation of center frequency specified! So I'm using ",self.glova.f0* 1e-12," THz")
            Df = 0
    
        self.Df = Df 
    
            
        if B == None and self.B != None:            
            B = self.B
                    
        if B == None and self.B == None:            
            B = 100.0  
            print ("WARNING: pypho_optfi: No FWHM bandwidth specified! So I'm using ", B, " GHz")    
                    
        self.B = B
        
        
        if loss == None and self.loss != None:
            loss = self.loss
        if loss == None and self.loss == None:            
            loss = 0
            print ("WARNING: pypho_optfi: No filter loss specified! So I'm using ", loss, " dB")                
            
        self.loss = loss
            

########################################################################
