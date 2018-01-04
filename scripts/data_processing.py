# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 10:52:36 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""
import pandas as pd
import numpy as np

def load_and_preproces_data(inputdata,taxon,filter_parameters,variables,res,nan_value):
    """ Load data and variables list
    
    Parameters
    ----------  
            'inputdata' (str): name of inputdata 
            'variables' (str): name of variables
            'taxon' (str): name of the taxon
            'res' (str): name of directory to write output
    
    Returns
    -------
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value"]
        'variables' (pandas df): List of considered variables
                            columns: ["variable","name_sim","consider"]                            
    """
    "Load inputdata"
    [inputdata,variables,filter_parameters] = load_data([inputdata,variables,filter_parameters])
    variables = variables["name_sim"][variables["consider"]==1].unique().tolist()
    
    "Fix nan in inputdata"
    inputdata = fix_nan(inputdata,nan_value)
    
    "Extract data for considered taxon and variables"
    inputdata = extract(inputdata,variables,taxon)
    from resample import resample_data
    
    " In earlier versions of code, automatic resampling option was available"
    " However, now resampling should be done outside SDMIT"
    "Part of code is left to assign unique ID's to each biological sample"

    resample = "False"
    inputdata = resample_data(inputdata,resample)

    "Extract parameters for considered variables"
    filter_parameters = filter_parameters[filter_parameters["variable"].isin(variables)]
    
    return inputdata,filter_parameters,variables

def load_data(files):
    """Load multiple ".csv" files in pandas dataframes
    
    Parameters
    ---------- 
        'files' (list): Files (str) which should be read to pandas dataframes
    
    Returns
    -------
        'data' (list): Pandas dataframes                        
    """
    
    data = [0]*len(files)
    
    for i in range(len(files)):
        
        data[i] = pd.read_csv(files[i],encoding = "ISO-8859-1")
        
    return data

def extract(inputdata,variables,taxon):
    """Extract/filter inputdata for the specified variables and taxon
    
    Parameters
    ---------- 
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value"]

        'variables' (pandas df): List of considered variables
                            columns: ["variable","name_sim',"consider"] 
        'taxon' (str): name of the taxon    
        
    Returns
    -------
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value"]               
    """
    
    "Extract data for taxon and variables "
    inputdata = inputdata[inputdata["variable"].isin(variables)]
    inputdata = inputdata[inputdata["taxon"]==taxon]
    "Remove NaN values"
    inputdata = inputdata[~inputdata["value"].isnull()]
    #inputdata["value"] = inputdata["value"].astype(float)
    
      
    return inputdata

def fix_nan(inputdata,nan_value):
    """insert nan_value for empty records in X, Y and date
    
    Parameters
    ---------- 
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value"]

        'nan_value' (float): value for nan value    
        
    Returns
    -------
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value"]             
    """
    col = ["date","X","Y"]

    for i in col:
        
        if np.sum(inputdata[i].isnull())>0:
        
            inputdata.loc[:,i] = inputdata.loc[:,i].fillna(-nan_value)
       
    return inputdata
#
#def prepare_dynamic_grid(filter_parameters,inputdata,runs):
#    
#    un_var = filter_parameters["variable"].unique()
#    
#    filter_parameters["lower_b"] = 0.
#    filter_parameters["upper_b"] = 0.
#    for i in un_var:
#
#        b2 = filter_parameters.loc[filter_parameters["variable"]==i,"b2"].values[0]
#        data_i = inputdata["value"][(inputdata["variable"]==i)]
#        filter_parameters.loc[filter_parameters["variable"]==i,"lower_b"] = Grid([np.percentile(data_i[inputdata["value"]>b2],i) for i in np.arange(0,runs,1)])
#        filter_parameters.loc[filter_parameters["variable"]==i,"upper_b"] = Grid([np.percentile(data_i[inputdata["value"]<b2],i) for i in np.arange(0,runs,1)])
#    
#    return filter_parameters
#    
#class Grid():
#    
#    def __init__(self,values):
#        
#        self.values = values
#    
#    def sample(self,n):
#        
#        return self.values[n]
#        
