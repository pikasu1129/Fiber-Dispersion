# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 17:13:42 2018

@author: arne
"""
# cython speedfiber.pyx -a
# python3 speedup.py build_ext --inplace

#defining NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

from __future__ import division

cimport numpy as np
import numpy as np

#import pyfftw
import time
import sys
import os

import cython
from cython.parallel import prange
from cython.view cimport array as cvarray
from cython import boundscheck, wraparound

cdef extern from "math.h" nogil:
    double sin(double x)
    double cos(double x)
    double log(double x)
    double sqrt(double x)
    double exp(double x)
    double abs(float x)
    double pow(double x, double y)
    double copysign( double x, double y)

cdef extern from "/usr/local/include/fftw3.h":
    int fftw_init_threads()
    void fftw_plan_with_nthreads(int)
   
    cdef int FFTW_FORWARD
    cdef int FFTW_BACKWARD
    cdef unsigned FFTW_MEASURE
    cdef unsigned FFTW_ESTIMATE
    cdef unsigned FFTW_PATIENT
    cdef unsigned FFTW_DESTROY_INPUT
    cdef unsigned FFTW_UNALIGNED
 
    ctypedef double fftw_complex[2]

    void *fftw_malloc(size_t)
    void fftw_free(void *)

    ctypedef struct _fftw_plan:
       pass

    ctypedef _fftw_plan *fftw_plan

    void fftw_execute(fftw_plan)
    void fftw_destroy_plan(fftw_plan)
    fftw_plan fftw_plan_dft_1d(int, fftw_complex*, fftw_complex*, int, unsigned)
    #void fftw_print_plan(byte)
    int fftw_import_system_wisdom();
    #int fftw_import_wisdom_from_string(const char *input_string)
    #int fftw_import_wisdom_from_file(const char *filename);
    int fftw_import_wisdom_from_filename(const char *filename)
    char *fftw_export_wisdom_to_string()
    int fftw_export_wisdom_to_filename(const char *filename)
    #void fftw_forget_wisdom()

cdef extern from "complex.h" nogil:
    double cabs(double complex)    
    double complex cexp(double complex)
    double creal( double complex z );
    double cimag( double complex z );
    

@boundscheck(False)
@wraparound(False)
def cyfiber(int noofsamples, float l, complex[:] E_x, complex[:] E_y, float alpha, float gamma_intern, float phi_max,
            double[:,:] birefarray, int birefsteps, float max_step, double beta_2, double beta_3, double frange, double fres, int fftw3_threads, fftw3_wisdom_path,
            complex[:] Ex_out, complex[:] Ey_out):

    cdef double noofsample_float = 1.0
    noofsample_float /= noofsamples
    fftw_init_threads()
    fftw_plan_with_nthreads(fftw3_threads)
    cdef int n = E_x.shape[0]
    cdef fftw_complex *Elx
    cdef fftw_complex *Ely
    cdef fftw_complex *outx
    cdef fftw_complex *outy
    Elx  = <fftw_complex *> fftw_malloc(n * sizeof(fftw_complex))
    Ely  = <fftw_complex *> fftw_malloc(n * sizeof(fftw_complex))
    outx = <fftw_complex *> fftw_malloc(n * sizeof(fftw_complex))
    outy = <fftw_complex *> fftw_malloc(n * sizeof(fftw_complex))
    # Wisdom

    erg = fftw_import_wisdom_from_filename(bytes(fftw3_wisdom_path, 'utf-8'))

    cdef fftw_plan plan_yf = fftw_plan_dft_1d(n, <fftw_complex *>Ely, <fftw_complex *>outy, FFTW_FORWARD,  (FFTW_PATIENT | FFTW_DESTROY_INPUT | FFTW_UNALIGNED))
    cdef fftw_plan plan_xf = fftw_plan_dft_1d(n, <fftw_complex *>Elx, <fftw_complex *>outx, FFTW_FORWARD,  (FFTW_PATIENT | FFTW_DESTROY_INPUT | FFTW_UNALIGNED))
    cdef fftw_plan plan_xb = fftw_plan_dft_1d(n, <fftw_complex *>outx, <fftw_complex *>Elx, FFTW_BACKWARD, (FFTW_PATIENT | FFTW_DESTROY_INPUT | FFTW_UNALIGNED))
    cdef fftw_plan plan_yb = fftw_plan_dft_1d(n, <fftw_complex *>outy, <fftw_complex *>Ely, FFTW_BACKWARD, (FFTW_PATIENT | FFTW_DESTROY_INPUT | FFTW_UNALIGNED))

    fftw_export_wisdom_to_filename(bytes(fftw3_wisdom_path, 'utf-8'));

    cdef float lc = 0
    cdef float powmax_x = 0.0
    cdef float powmax_y
    cdef double[:] beta_xr = np.zeros(noofsamples)   
    cdef double[:] beta_yr = np.zeros(noofsamples) 
    cdef double[:] beta_xi = np.zeros(noofsamples)
    cdef double[:] beta_yi = np.zeros(noofsamples) 

    cdef complex beta_x  
    cdef complex beta_y
    
    cdef double[:] power_x = np.zeros(noofsamples)
    cdef double[:] power_y = np.zeros(noofsamples)

    
    cdef double nl_x 
    cdef double nl_y   
    cdef double Exr_tmp, Exi_tmp, Eyr_tmp, Eyi_tmp, Ex_tmp, Ey_tmp
    cdef double Exr_tmp2, Eyr_tmp2
    
    cdef double power1_x
    cdef double power1_y

    cdef double[:] E_x_r = np.zeros(noofsamples)
    cdef double[:] E_x_i = np.zeros(noofsamples)
    cdef double[:] E_y_r = np.zeros(noofsamples)
    cdef double[:] E_y_i = np.zeros(noofsamples)

    
    cdef double E_tmp

    cdef double cos_x
    cdef double sin_x
    cdef double cos_y
    cdef double sin_y
    
    #cdef complex[noofsamples] beta_term
    cdef complex delta_beta
    cdef complex delta_beta_current

    cdef double[::1] zleff      = -np.log( np.array([0.98, 0.95, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.4, 0.3, 0.25, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 1e-10]) ) / (alpha + 1.0e-12) #kleiner pfusch
    
    cdef int zleff_i            = 0
    cdef doCoordRot             = False
    cdef doCalcLeff             = True
    cdef float next_l           = 0
    cdef float next_leff        = 0
    cdef float fak_xpm          = 2.0/3.0

    cdef float gamma_fast       = -0.0000000001
    cdef int birefindex         = -1
    cdef double cos_angle       = 0.0
    cdef double sin_angle       = 0.0
    cdef double biref_ang       = 0.0
    cdef double  Domega_max      = frange * 1.0 * 3.141592653589793  #/ 1.0e12
    cdef double  Domega_diff     = 2.0 * 3.141592653589793 * fres # / 1.0e12
    cdef double  Domega, beta_fac, 
    cdef complex var_tmp_1, var_tmp_2    
    


    cdef int i

    for i in prange(noofsamples, nogil=True):
        Elx[i][0] = creal(E_x[i])
        Elx[i][1] = cimag(E_x[i])
        Ely[i][0] = creal(E_y[i])
        Ely[i][1] = cimag(E_y[i])
    
    del(E_x, E_y)
    print('Start ciber...', alpha)
    tic0 = time.time()
    while lc < l :
        
        if birefindex < birefsteps and lc + next_l >= birefarray[birefindex][1]:
            #print( lc + next_l , birefindex, birefarray[birefindex].z_point)
            doCoordRot = True
            
        # TODO: stimmt das?
        if (lc + next_l) >= l:
             doCalcLeff = True

        if lc >= zleff[zleff_i] or doCoordRot or doCalcLeff:               # calculate the steplength
            doCalcLeff = False
            #print (lc, zleff[zleff_i], zleff_i, sqrt(exp(-alpha * next_l)))
            if lc >= zleff[zleff_i] :
                zleff_i +=1
            
            powmean_x = 0
            powmean_y = 0
            for i in range(noofsamples):        # das kann man noch verbessern!
                powmean_x += Elx[i][0]*Elx[i][0] + Elx[i][1]*Elx[i][1]
                powmean_y += Ely[i][0]*Ely[i][0] + Ely[i][1]*Ely[i][1]            
            powmean_x /= noofsamples
            powmean_y /= noofsamples
            powmax_x = np.mean([powmean_x, powmean_y])
            next_leff = phi_max / powmax_x
            
            if next_leff*alpha < 1.0 and alpha > 0.0:
                next_l = -log(1.0 - next_leff*alpha)/alpha
            else:
                next_l = next_leff

            if next_l > max_step:
                next_l = max_step
                if alpha > 0.0:
                    next_leff = (1.0 - exp(-alpha*next_l)) / alpha
                elif alpha == 0.0:
                    next_leff = next_l

            if (lc + next_l) >= l:
                next_l =  l - lc
                if alpha > 0.0:
                    next_leff = (1.0 - exp(-alpha*next_l)) / alpha
                elif alpha == 0.0:
                    next_leff = next_l

            doCoordRot = False
            if birefindex < birefsteps and lc + next_l >= birefarray[birefindex][1]:
                if birefarray[birefindex][1] > 0 :
                    next_l = birefarray[birefindex][1] - lc
                    if alpha > 0.0:
                        next_leff = (1.0 - exp(-alpha*next_l)) / alpha
                    elif alpha == 0.0:
                        next_leff = next_l
                delta_beta_current = birefarray[birefindex][2]
                birefindex = birefindex + 1
                doCoordRot = True
                doCalcLeff = True

            var_tmp_1 = 0.5*delta_beta_current * next_l * 1.0e-12
            var_tmp_2 = sqrt(exp(-alpha * next_l))
           

            #for i in range(noofsamples):
            for i in prange(noofsamples, nogil=True):  
                Domega =  i*Domega_diff - (copysign(1., i-noofsamples/2) + 1.0) * Domega_max
                beta_fac   = beta_2*0.5 * pow(Domega, 2.0) + beta_3 / 6.0  * pow(Domega, 3.0) 
                
                beta_x = cexp(-1.0j * beta_fac * next_l) * var_tmp_2 * noofsample_float# for x- and y-pol
                beta_y = beta_x * cexp( +1.0j*Domega * var_tmp_1)# x-pol # hier nur 1.0 auf eine Polarisation. Selber Effekt, aber weniger variablen.
                beta_x = beta_x * cexp( -1.0j*Domega * var_tmp_1)# y-pol
                
                beta_xr[i] = creal(beta_x)
                beta_xi[i] = cimag(beta_x)

                beta_yr[i] = creal(beta_y)
                beta_yi[i] = cimag(beta_y)

            gamma_fast = -next_leff*gamma_intern
       
        # Nonlinear fiber effects
        tic = time.time() 
        for i in prange(noofsamples, nogil=True):
            
            Ex_tmp = Elx[i][0]*Elx[i][0] + Elx[i][1]*Elx[i][1]            
            Ey_tmp = Ely[i][0]*Ely[i][0] + Ely[i][1]*Ely[i][1]
            
            nl_x = (Ex_tmp  + Ey_tmp * fak_xpm) * gamma_fast                   
            nl_y = (Ey_tmp  + Ex_tmp * fak_xpm) * gamma_fast
            
            cos_x = cos(nl_x) #1.0 - nl_x*nl_x * (0.5 - nl_x*nl_x * 0.041666666666666664)
            sin_x = sin(nl_x) #nl_x - nl_x*nl_x*nl_x * (0.16666666666666666 - nl_x*nl_x * 0.008333333333333333)
            cos_y = cos(nl_y) #1.0 - nl_y*nl_y * (0.5 - nl_y*nl_y * 0.041666666666666664)
            sin_y = sin(nl_y) #nl_y - nl_y*nl_y*nl_y * (0.16666666666666666 - nl_y*nl_y * 0.008333333333333333)

            Ex_tmp = Elx[i][0]
            Elx[i][0] = Ex_tmp * cos_x - Elx[i][1] * sin_x
            Elx[i][1] = Ex_tmp * sin_x + Elx[i][1] * cos_x
            Ey_tmp = Ely[i][0]
            Ely[i][0] = Ey_tmp * cos_y - Ely[i][1] * sin_y
            Ely[i][1] = Ey_tmp * sin_y + Ely[i][1] * cos_y

        
        
        # Linear fiber effects
        fftw_execute(plan_xf)   # FFT fw x-pol
        fftw_execute(plan_yf)   # FFT fw y-pol
        
        for i in prange(noofsamples, nogil=True):

            Exr_tmp = outx[i][0]
            outx[i][0] = Exr_tmp*beta_xr[i] - outx[i][1]*beta_xi[i] # x-pol linstep real                       
            outx[i][1] = outx[i][1]* beta_xr[i] + Exr_tmp*beta_xi[i] # x-pol linstep imag
            
            Eyr_tmp = outy[i][0]
            outy[i][0] = Eyr_tmp*beta_yr[i] - outy[i][1]*beta_yi[i] # y-pol linstep real
            outy[i][1] = outy[i][1]* beta_yr[i] + Eyr_tmp*beta_yi[i] # y-pol linstep imag
        
        fftw_execute(plan_xb)  # FFT bw x-pol
        fftw_execute(plan_yb)  # FFT bw y-pol

        # Polarisation
        if doCoordRot:
            biref_ang += birefarray[birefindex][0]
            sin_angle = sin(birefarray[birefindex][0])
            cos_angle = cos(birefarray[birefindex][0])
            for i in prange(noofsamples, nogil=True):
                Exr_tmp = Elx[i][0]
                Exi_tmp = Elx[i][1]
                Eyr_tmp = Ely[i][0]
                Eyi_tmp = Ely[i][1]
                Elx[i][0] = Exr_tmp * cos_angle - Eyr_tmp * sin_angle
                Elx[i][1] = Exi_tmp * cos_angle - Eyi_tmp * sin_angle
                Ely[i][0] = Exr_tmp * sin_angle + Eyr_tmp * cos_angle
                Ely[i][1] = Exi_tmp * sin_angle + Eyi_tmp * cos_angle
            doCoordRot = False

        lc += next_l
        curti = time.time()
        sys.stdout.write( "\r[%s%s] %5.2f%%, %8.2f m, %5.1f sec, %5.2f m, %6.1f steps/sec" % ( '#'*int(50*lc/l), '.'*(50-int(50*lc/l)), 100*lc/l, lc ,(curti - tic0 ) , next_l, 1.0/(curti - tic)  ) )
        
    fftw_destroy_plan(plan_xf)
    fftw_destroy_plan(plan_xb)
    fftw_destroy_plan(plan_yf)
    fftw_destroy_plan(plan_yb)
    
    sin_angle = sin( -biref_ang )
    cos_angle = cos( -biref_ang )
    for i in prange(noofsamples, nogil=True):
        Exr_tmp = Elx[i][0]
        Exi_tmp = Elx[i][1]
        Eyr_tmp = Ely[i][0]
        Eyi_tmp = Ely[i][1]
        Elx[i][0] = Exr_tmp * cos_angle - Eyr_tmp * sin_angle
        Elx[i][1] = Exi_tmp * cos_angle - Eyi_tmp * sin_angle
        Ely[i][0] = Exr_tmp * sin_angle + Eyr_tmp * cos_angle
        Ely[i][1] = Exi_tmp * sin_angle + Eyi_tmp * cos_angle
    
    #print('Ende loop!')
   
    for i in prange(noofsamples, nogil=True):
        Ex_out[i] = Elx[i][0] + 1j*Elx[i][1]
        Ey_out[i] = Ely[i][0] + 1j*Ely[i][1]

    fftw_free(Elx)
    fftw_free(Ely)
    fftw_free(outx)
    fftw_free(outy)
    # TODO : Destroy plans fftw_destroy_plan(plan);
    #print('Ende ciber!')
    
    
