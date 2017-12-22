# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 09:22:33 2014
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import time
import datetime
import numpy as np

class Runtime():
    
    def __init__(self,n,delta_percentage):

        # initialize percentages to evaluate 
        self.p = range(5,100+int(delta_percentage),int(delta_percentage))
        # initialize time matrix
        self.time = [0.]*(int(n)+1);
        # put in first timing
        self.time[0] = time.time();self.t0=time.time()
        # hold iterator for self.p
        self.iter = 0;
        
    def iteration(self,i):

        # calculate current percentrage
        current_percentage=float(i)/float(len(self.time))*100.

        # if perc is over the fixed one
        if self.p[self.iter]<current_percentage:
            # time to do
            dt=datetime.timedelta(seconds=np.nanmean(self.time[1:i])*(len(self.time)-i))
            print(str(int(current_percentage))+"% done... "+str(dt)+" remaining ...")
            # increase iterator for self.p            
            self.iter+=1  
        else:
            # calculate difference in time
            self.time[i] = time.time()-self.time[0];
            # save current time in first record of self.time
            self.time[0] = time.time()
             
    def close(self,print_value=True):
        
        dt=datetime.timedelta(seconds=time.time()-self.t0)

        if print_value==True:

            print("End simulation, total runtime is "+str(dt))
        
        return dt
        
