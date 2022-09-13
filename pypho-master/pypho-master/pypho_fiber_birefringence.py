# -*- coding: utf-8 -*-

########################################################################

import numpy as np

"""
Class for modeling a birefringent fiber. To model birefringent behavior, the 
pypho_fiber class must be given an array of this class as an argument for
"birefarray".

Each pypho_fiber_birefringence object represents a segment of the fiber. The 
z_point variable is the z-coordinate of the start of the segment. The first
z-coordinate in the array must always be 0.

You also have to provide the angle by which the polarisation axis is rotated
and the delta beta for the birefringence.
"""

class pypho_fiber_birefringence(object):
    def __init__(self, z_point = None, angle = None, delta_beta = None):
        
        if z_point == None:
            self.z_point = 0.0
        else:
            self.z_point = z_point
            
        if angle == None:
            self.angle = 0.0
        else:
            self.angle = angle
            
        if delta_beta == None:
            self.delta_beta = 0.0
        else:
            self.delta_beta = delta_beta
            
    """Method for automatically creating a pypho_fibre_birefringence array with
    randomized delta beta and rotation angles."""
    @staticmethod
    def create_pmd_fibre (totalLength, stepLength, maxDeltaBeta):
        currentlength = 0.0
        fibres = []
        while currentlength < totalLength:
            fibres.append(pypho_fiber_birefringence(currentlength, np.random.rand()*2*np.pi, np.random.rand()*maxDeltaBeta))
            currentlength += stepLength
            
        return fibres