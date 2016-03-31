# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:13:14 2016

@author: sacha
"""

import pandas as pd
pd.set_option('chained_assignment',None)
import numpy as np
import os

def resampling_protocol(inputdata,fields,res,full_output=False):

    """ 
    Resampling protocol to sample for model development and optimisation
    NOTE: Data are resampled so the P/A-ratio, class-ratio (biological index),
            and waterbody type are equal!
    Arguments:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value"]

        'res' (str): name of directory to write output  
        
        'full_output' (bool): self-explanatory
        
    Returns:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value","development"]                     
    """
    
    inputdata = algorithm(inputdata,fields)
    if full_output==True:
        write(inputdata,fields,res) 
    
    return inputdata
    
def algorithm(inputdata,fields):
    """
    Algorithm to sample for model development and optimisation$
    Arguments:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","X","Y",[fields],"taxon","abundance","variable","value"]

        'fields' (list): Names of columns on which the sampling is done. The data are sampled 
                         so all classes in the fields are equally distributed over the samples
        
    Returns:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","X","Y",[fields],"taxon","abundance","variable","value","development"]      
    """
    
    " Initiate IDs dictionary to save IDs for model development and optimisation"
    folds = ["development","optimisation"]
    inputdata["fold"] = np.nan
    
    IDs = {i:[] for i in  folds}

    conds = ['inputdata["abundance"]!=0','inputdata["abundance"]==0']    
    
    for i in conds:
        
        sample = inputdata[["ID","X","Y","taxon","abundance"]+fields][eval(i)].drop_duplicates()
        "Sample in space according to the defined fields"
        IDs = spatial_sampler(IDs,sample,fields)
    
    "Appoint whether points are used for model development or optimisation "
    for i in folds:

        inputdata["fold"][inputdata["ID"].isin(IDs[i])] = i
    
    return inputdata
    
def spatial_sampler(IDs,sample,fields):
    """
    Sample IDs according to the defined fields and the X and Y coordinate
    Arguments:
    
        'IDs' (dictionary): lists of ID instances per fold (keys)
        'sample' (pandas df): sample of Biological and environmental measurements
                            columns: ["ID","X","Y",[fields],"taxon","abundance","variable","value"]
        
    Returns:
    
        'IDs' (dictionary): updated lists of ID instances per fold (keys)
    """
    
    "Sample first point "
    point = sample.sample()
    "Get distance "
    sample["d"] = np.sqrt((sample["X"]-point["X"].values[0])**2+(sample["Y"]-point["Y"].values[0])**2)
    "Sort according to fields and distance"
    sample =  sample.sort(fields+["d"]).reset_index()
    
    for i in range(len(IDs.keys())):
        
        "Sample each alternate point"
        r = range(i,len(sample),2)
        IDs[IDs.keys()[i]] = IDs[IDs.keys()[i]]+sample["ID"].iloc[r].tolist()

    return IDs
    
def write(inputdata,fields,res):
    
    "Save folds" 
    inputdata.to_csv(os.path.join(res,"sample_data.csv"))
    
    "Print distributions of samples"
    inputdata["number of samples"] = 1
    inputdata["abundance"][inputdata["abundance"]!=0] = 1
    inputdata = inputdata[["ID","fold","abundance"]+fields+["number of samples"]].drop_duplicates()
    inputdata.groupby(["fold","abundance"]+fields).aggregate({"number of samples":np.sum}).reset_index().to_csv(os.path.join(res,"distribution_samples.csv"))
    
