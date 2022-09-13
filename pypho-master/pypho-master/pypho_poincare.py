# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 11:14:36 2018

@author: optsim1
"""

# Creating Poincare Sphere

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pypho_setup import pypho_setup
from matplotlib.legend_handler import HandlerLine2D
import sys



########################################################################
class pypho_poincare(object):
    
    def __init__(self, glova = None, color = None, marker = None):
         
        if glova == None:
            print ("ERROR: pypho_optfi: You must define the global variables")
            sys.exit("PyPho stopped!")
        self.glova  = glova
        self.E      = None
        self.color  = None
        self.marker = None
        
        self.set(self.E, color, marker)
 
########################################################################       
    def __call__(self, E = None, color = None, marker = None):
        
        if len(E) <> 2:
            print ("ERROR: You must define an optical signal")
            sys.exit("PyPho stopped!")
            
        self.set(E, color, marker)
        self.out(E)
 
########################################################################       
    def out(self, E = None):
        
        s0_param, s1_param, s2_param, s3_param = self.calc_stokes_param(E[0,:], E[1, :])                                                                                            
        
        fig = plt.gcf()
        poincare_plt = fig.add_subplot(111, projection='3d')
        poincare_plt.set_aspect("equal")
        
        poincare_plt.set_xlabel('$S_2$', fontsize = 20)
        poincare_plt.set_ylabel('$S_1$', fontsize = 20)
        poincare_plt.set_zlabel('$S_3$', fontsize = 20, rotation = 0)
        poincare_plt.dist = 6


##################Wireframe Plot#######################################################################################      

        r = 1
        pi = np.pi
        cos = np.cos
        sin = np.sin
        phi, theta = np.mgrid[0.0:pi:50j, 0.0:2.0*pi:50j]
        x = r*sin(phi)*cos(theta)
        y = r*sin(phi)*sin(theta)
        z = r*cos(phi)        
        
        poincare_plt.plot_surface( x, y, z,  rstride=1, cstride=1, color='c', alpha=0.1, linewidth=0)        
        
        w = 1.0
        poincare_plt.plot([0, 0],[-w, w], [0, 0], '0.3', linewidth = 0.25)
        poincare_plt.plot([-w, w],[0, 0], [0, 0], '0.3', linewidth = 0.25)
        poincare_plt.plot([0, 0],[0, 0], [-w, w], '0.3', linewidth = 0.25)
        
        w = 0.2
        poincare_plt.plot([0, 0],[0, -w], [0, 0], '0.0', linewidth = 0.5)
        poincare_plt.plot([0, w],[0, 0], [0, 0], '0.0', linewidth = 0.5)
        poincare_plt.plot([0, 0],[0, 0], [0, w], '0.0', linewidth = 0.5)
        
        w=0.25
        poincare_plt.text(w, 0, 0, "$S_2$", color='0.35')
        poincare_plt.text(0, -w, -0.1, "$S_1$", color='0.35')
        poincare_plt.text(0, 0, w, "$S_3$", color='0.35')    
        
        # Cylinder
        x=np.linspace(-1, 1, 100)
        z=np.linspace(-0.01, 0.01, 10)
        Xc, Zc=np.meshgrid(x, z)
        Yc = np.sqrt(1-Xc**2)
        
        # Draw parameters
        poincare_plt.plot_surface(Xc,  Yc, Zc, alpha=0.5, color='c')
        poincare_plt.plot_surface(Xc, -Yc, Zc, alpha=0.5, color='c')

        poincare_plt.plot_surface( Zc, Xc, Yc, alpha=0.5, color='c')
        poincare_plt.plot_surface( Zc, Xc,-Yc, alpha=0.5, color='c')

        poincare_plt.plot_surface( Xc, Zc, Yc, alpha=0.5, color='c')
        poincare_plt.plot_surface( Xc, Zc,-Yc, alpha=0.5, color='c')        
        
        col = '0'
        w=1.1
        poincare_plt.text(0, 0, +w, "$R\circlearrowright$", color=col, bbox={'facecolor':'white', 'alpha':0.5, 'pad':2})    
        poincare_plt.text(0, 0, -w, "$L\circlearrowleft$", color=col, bbox={'facecolor':'white', 'alpha':0.5, 'pad':2})    
        poincare_plt.text(0, -w, 0, "$H\leftrightarrow$", color=col, bbox={'facecolor':'white', 'alpha':0.5, 'pad':2})    
        poincare_plt.text(0, +w, 0, "$V\updownarrow$", color=col, bbox={'facecolor':'white', 'alpha':0.5, 'pad':2})    
        poincare_plt.text(+w, 0, 0, "$D\\nearrow$", color=col, bbox={'facecolor':'white', 'alpha':0.5, 'pad':2})    
        poincare_plt.text(+w, 0, 0, "$D\swarrow$", color=col, bbox={'facecolor':'white', 'alpha':0.5, 'pad':2})   
        poincare_plt.text(-w, 0, 0, "$A\\nwarrow$", color=col, bbox={'facecolor':'white', 'alpha':0.5, 'pad':2})    
        poincare_plt.text(-w, 0, 0, "$A\searrow$", color=col, bbox={'facecolor':'white', 'alpha':0.5, 'pad':2})           
        
        poincare_plt.axis("off")
     
        

  
##################Scatter Plot####################################################################################### 

        for xp in range(np.size(s0_param)): 
           poincare_plt.scatter([s2_param[xp]], [-s1_param[xp]], [s3_param[xp]], color = self.color ,alpha = 1, marker = self.marker, s = 50)
           xp += 1
        
        
        fig.show()
        
        



########################################################################        
    def calc_stokes_param(self, A_x, A_y):
        
        S0_param = (np.abs(A_x))**2 + (np.abs(A_y))**2
        S1_param = ((np.abs(A_x))**2 - (np.abs(A_y))**2)/S0_param
        S2_param = (2*np.abs(A_x)*np.abs(A_y)*np.cos(np.angle(A_x)-np.angle(A_y)))/S0_param
        S3_param = (2*np.abs(A_x)*np.abs(A_y)*np.sin(np.angle(A_x)-np.angle(A_y)))/S0_param    
        
        return(S0_param, S1_param, S2_param, S3_param)
        
 
########################################################################
    def set(self, E = [], color = None, marker = None):


        if color == None and self.color == None:
             self.color = "r"
        elif color != None:
            self.color = color

        if marker == None and self.marker == None:
             self.marker = "."
        elif marker != None:
            self.marker = marker

########################################################################      
    
    