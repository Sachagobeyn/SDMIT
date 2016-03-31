# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 10:52:36 2015

@author: sacha
"""

import pandas as pd
pd.set_option('chained_assignment',None)

def load_and_sample_data(inputdata,taxon,variables,res):
    """ 
    Load data and variables list and sample the data
    
    Arguments:
        
            'inputdata' (str): name of inputdata 
            NOTE: structure inputdata should be conform Section XX tutorial
            'variables' (str): name of variables
            NOTE: structure variables should be conform Section XX tutorial
            'taxon' (str): name of the taxon
            'res' (str): name of directory to write output
    
    Returns:
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value","development"]
        'variables' (pandas df): List of considered variables
                            columns: ["variable","consider"]                            
    """
    "Load inputdata"
    [inputdata,variables] = load_data([inputdata,variables])
    variables = variables["variable"][variables["consider"]==1].unique().tolist()
        
    "Extract data for considered taxon and variables"
    inputdata = extract(inputdata,variables,taxon)
    
    "Split data (if no development tag is defined in the inputdata)"
    if "fold" not in inputdata:
        
        from resampling_protocol import resampling_protocol
        return resampling_protocol(inputdata,["type","EQRclass"],res,full_output=True),variables
    
    else:
        
        return inputdata,variables
    
def load_data(files):
    """ 
    Load multiple ".csv" files in pandas dataframes
    
    Arguments:
        
        'files' (list): Files (str) which should be read to pandas dataframes
    
    Returns:
    
        'data' (list): Pandas dataframes                        
    """
    
    data = [0]*len(files)
    
    for i in range(len(files)):
        
        data[i] = pd.read_csv(files[i])
        
    return data

def extract(inputdata,variables,taxon):
    """ 
    Extract/filter inputdata for the specified variables and taxon
    
    Arguments:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value",optional="development"]

        'variables' (pandas df): List of considered variables
                            columns: ["variable","consider"] 
        'taxon' (str): name of the taxon    
        
    Returns:
    
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value",optional="development"]                     
    """
    
    "Extract data for taxon and variables "
    inputdata = inputdata[inputdata["variable"].isin(variables)]
    inputdata = inputdata[inputdata["taxon"]==taxon]
    "Remove NaN values"
    inputdata = inputdata[~inputdata["value"].isnull()]
    inputdata["value"] = inputdata["value"].astype(float)
    
      
    return inputdata
    
