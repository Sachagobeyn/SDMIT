# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:13:14 2016
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import pandas as pd
import random
from copy import deepcopy
import numpy as np

"""NOTE: main functionalities in function are not used as responsability of resampling
is r
"""
def resample_data(inputdata,resample=False):

    """ 
    Resampling protocol to sample for model development and optimisation
    NOTE: Data are resampled so the P/A-ratio, class-ratio (biological index),
            and waterbody type are equal!
    Arguments:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value"]

        'resample' (boolean): Yes (True) or No (False)
        
    Returns:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value"]                     
    """
    
    inputdata = algorithm(inputdata,resample)
    
    return inputdata
    
def algorithm(inputdata,resample):
    """
    
    NOTE: ONLY USED TO ASSIGN IDs
    Arguments:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","X","Y",[fields],"taxon","abundance","variable","value"]

        'fields' (list): Names of columns on which the sampling is done. The data are sampled 
                         so all classes in the fields are equally distributed over the samples
        
    Returns:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","X","Y",[fields],"taxon","abundance","variable","value"]      
    """
    
    "IDs"
    IDs = inputdata["ID"].unique().tolist()

    "Sample IDs with replacement"
    sample_IDs = [random.choice(IDs) for i in range(len(IDs))]
    
    "Select samples from data and return"
    counter = 0
    
    "Resample: yes/no"
    if eval(resample):
    
        import sys        
        print("[PROGRAMMED EXIT] \n\t"+"Resampling is not implemented in SDMIT version 2!")
        sys.exit("="*19)         

        for i in sample_IDs:

            s = inputdata[inputdata["ID"] == i]
            s["sample"] = counter
            if counter==0:
                sample = s
            else:
                sample = sample.append(s)
            counter += 1
            
    else:
        
        sample = deepcopy(inputdata)
        
        if "sample" not in inputdata:
            # couple sample number
            sampleID = pd.DataFrame(data=np.zeros([len(IDs),2]),columns=["ID","sample"])
            sampleID["sample"] = np.arange(len(IDs))
            sampleID["ID"] = IDs
            sample  = sample.merge(sampleID,how="left",on="ID")
        
    return sample
    
