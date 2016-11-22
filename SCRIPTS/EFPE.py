# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 11:03:22 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""
import os
import numpy as np
import pandas as pd
pd.set_option('chained_assignment',None)

def EFPE(inputdata,model_settings,taxon,res):
    """ 
    Environemntal filter parameter estimation
    
    Arguments:
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value","development"]      
        'model_settings' (dictionary): considered percentiles of the environmental 
                                    data to estimate parameters of the habitat
                                    preference curve.
                            
        'taxon' (str): name of taxon
        'res' (str): name of map to write output
    
    Returns:
        'parameters' (pandas df): Estimated parameters for habitat preference curves
                            columns: ["taxon","a1","a2","a3","a4","type"]
    """
    
    "Take data for development"
    sample = inputdata[inputdata["fold"]=="development"]
    "Take all presences into account"    
    sample = sample[sample["abundance"]!=0]
    "Estimate the parameters of the habitat preference curve"
    parameters = curve_parameter_estimation(sample,model_settings,taxon)  
    "Write"
    parameters.to_csv(os.path.join(res,"parameters_"+taxon+".csv"))
    
    return parameters
    
def curve_parameter_estimation(inputdata,model_settings,taxon):

    """ 
    Function to estimate habitat preference curves
    
    Arguments:
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value","development"]      
        'model_settings' (dictionary): considered percentiles of the environmental 
                                    data to estimate parameters of the habitat
                                    preference curve.
                            
        'taxon' (str): name of taxon
    
    Returns:
        'parameters' (pandas df): Estimated parameters for habitat preference curves
                            columns: ["taxon","a1","a2","a3","a4","type"]
    """
    
    for i in model_settings.keys():     
        inputdata[i] = inputdata["value"][:]
        
    cond = "{"+",".join(["'"+i+"':lambda x:np.percentile(x,"+str(model_settings[i])+")" for i in model_settings.keys()])+"}"
    parameters = inputdata.groupby(["taxon","variable"]).aggregate(eval(cond)).reset_index()
    parameters["type"] = "acute"
    parameters["taxon"] = taxon
    
    return parameters
    
#def curve_parameter_estimation_abundance(data,perc,taxa):
#    
#    a = ["a1","a2","a3","a4"]
#
#    un_var = data["variable"].unique().tolist()
#    col = a + ["variable","taxa"]
#    
#    p = pd.DataFrame(np.zeros([len(un_var),len(col)]),columns = col)
#    p["taxa"] = taxa
#    p["variable"] = un_var
#    
#    for i in un_var:
#        
#        cond = p["variable"] == i
#        # sort 
#        x = data[data["variable"]==i]
#        x = x.sort(["value"]).reset_index()
#        x["abundance"] = np.log(x["abundance"])
#        x["p"] = [np.sum(x["abundance"].iloc[:k+1])/np.sum(x["abundance"].iloc[:]) for k in range(0,len(x["abundance"]))]
#    
#        # interpolate
#        xp = x["value"] 
#        yp = x["p"]
#        
#        for j in range(len(perc)):
#            
#            p[a[j]][cond] = np.interp(perc[j], yp,xp)
#  
#    
#	p[a[0]][cond] = np.min(x["value"])
#    	p[a[3]][cond] = np.max(x["value"])
#      
#    p["type"] = "acute"
#    
#    return p
#    
