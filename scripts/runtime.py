# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 09:22:33 2014
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import time
import datetime

class Runtime():
    """ Class to calculate runtime
    
    Attributes
    -----------
        't0' (time object): start iteration
    """
    def __init__(self):

        """ initialisation
    
        Parameters
        ----------
            none
        
        Return
        ------
            none
        """
    
        # put in first timing
        self.t0=time.time()
        # hold iterator for self.p
        
             
    def close(self,print_value=True):
        """ close timer
        
        Parameters
        ----------
            print_value (boolean): print to screen
        
        Return
        ------
            none
        """        
        dt=datetime.timedelta(seconds=time.time()-self.t0)

        if print_value==True:

            print("End simulation, total runtime is "+str(dt))
        
        return dt
        
