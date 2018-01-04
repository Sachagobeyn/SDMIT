# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:39:36 2017

@author: sacha
"""

from scipy.optimize import fsolve
import numpy as np

def hyperparameters(l,FE):

    """ Approximate optimal hyper parameters for genetic algorithm
    Based on formula Gibbs, M.S., Dandy, G.C., Maier, H.R., 2008. A genetic 
    algorithm calibration method based on convergence due to genetic drift. 
    Inf. Sci. (Ny). 178, 2857–2869.
    
    Parameters
    ----------        
    'l' (float): (maximum) size of genotype (equal to 4 multiplied by number of 
    variables (*for embedded feature selection!!!*) OR equal to the number of 
    variables (*for wrapper feature selection*)).
    
    'FE' (float):  Function evaluation, determined by dividing the computer time 
    available by the average time to compute the fitness function.
        
    Returns
    -------
    'N' (float): an approximate of the population size
    'm' (float): an approximate of the mutation rate
    'FE/N' (float): an approximate of the (minimum) number of iteration cycles
    """

    M = 3
    
    f = lambda x:FE/x*np.log10(1.-1./x)+M+np.log10(np.sqrt(l/12))
    
    N = fsolve(f,10)[0]
    N = -N if N<0 else N
    pm = 5/N
    
    return {'PS':N,'pm':pm,'iterations':FE/N}
    
