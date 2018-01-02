# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:17:52 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import numpy as np
from copy import deepcopy
#import sys

# x = predicitions
# y = observations

def calculate_performance(model_output,K,w=0.5,evaluate=False,Kmax=np.nan,threshold=np.nan):
    """ 
    Calculate peformance based on the simulated HSI and observed presence/absence (abundance)
    NOTE: Multiple self-defined formula's can be defined and saved in the dictionary 'criteria'
    
    Arguments:
        'model_output' (pandas df), columns
            'HSI' (float): Simulations
            'abundance'(float): Observed abundances                
        'variables' (list): self-explanatory
        
    
    Returns:
        'criteria' (dictionary): Dictionary with self-defined criteria
        NOTE: One of the criteria defined here has to be defined in ga_settingsfile!
 
    """
    
    criteria = {}
    "extract the predictions (x) and observations (y)"
    x = model_output["prediction"].values
    y = model_output["observation"].values
    
    "SSE"
    criteria["SSE"] = SSE(x,y)
    #criteria["wSSE"] = wSSE(x,y,w)
    
    "AIC (corrected for small sample size"
    N = len(model_output)
    _,criteria["AIC"] = nAIC(criteria["SSE"],N,K)
    criteria["BIC"] = nBIC(criteria["SSE"],N,K)
    #_,criteria["AIC"] = AIC(criteria["SSE"]/N,N,K)
    
    criteria["K"] = K
    criteria["N"] = N
    
    if evaluate==True:
        
        "Calculate Kappa, CCI, Sn, Sp and TSS (choose threshold value that maximises Kappa)"           
        performance = performance_threshold(x,y,["Kappa","Sn","Sp","TSS","CCI"],threshold=threshold,evaluation_criterion="TSS")
        criteria["Kappa"] = performance["Kappa"] 
        criteria["Sn"] = performance["Sn"] 
        criteria["Sp"] = performance["Sp"] 
        criteria["TSS"] = performance["TSS"]  
        criteria["CCI"] = performance["CCI"] 
        criteria["threshold"] = performance["threshold"]
        
        "Calculate AUC"
        from sklearn.metrics import roc_curve,auc
        fpr, tpr, _ = roc_curve(y,x)
        criteria["AUC"] = auc(fpr,tpr)

    return criteria 
    
def performance_threshold(x,y,criteria,threshold=np.nan,evaluation_criterion="TSS"):
    """Note: we really need to clean this code up"""
    performance = {i:0. for i in criteria}
    
    "No threshold defined: Get best from one of the defined criteria"
    if type(threshold) is not str:
        
        if np.isnan(threshold):
            
            l = np.arange(0,1.05,0.05)
            v = np.zeros([len(l)])
    
            "Check which threshold is the best, given the evalution criterion"        
            for i in range(len(l)):
                
                x_i = deepcopy(x)
                x_i[x_i<l[i]] = 0
                x_i[x_i!=0] = 1
                v[i] = eval(evaluation_criterion+"(x_i,y)")
                
            opt_threshold = np.nanmean(l[v==np.nanmax(v)])
        
        else:
            
            opt_threshold = threshold
    
    else:
        opt_threshold = threshold
    
    performance["threshold"] = opt_threshold
    
    "Get results with the found or defined threshold"
    for i in criteria:
        
        x_i = deepcopy(x)    
        
        if type(threshold) is str:
            
            if threshold=="prob":
    
                x_i = np.random.binomial(1, x_i)
        else:
            x_i[x_i<=opt_threshold] = 0
            x_i[x_i!=0] = 1
            
        performance[i] = eval(i+"(x_i,y)")        
    
    return performance
    
         
def Kappa(X,Y):
    
    a = np.sum(X+Y==2)
    b = np.sum(X-Y==1)
    c = np.sum(Y-X==1)
    d = np.sum(X+Y==0)
    
    tot = a + b + c + d
    Pa  = float(a + d)/tot
    PA1 = float(a + b)/tot
    PA2 = 1.0- PA1
    PB1 = float(a + c) /tot
    PB2 = 1.0 -PB1
    Pe  = PA1 *PB1 + PA2*PB2
    
    return (Pa -Pe)/ (1.0 -Pe)
    
def CCI(X,Y):
    
    return np.sum(X==Y)/float(len(X))
    
# sensitivity (Sn): correclty classified as present
def Sn(x,y):
    
    if float(sum(x+y==2) + sum(x-y==-1))!=0:
    
        return float(sum(x+y==2)) / float(sum(x+y==2) + sum(x-y == -1))
    
    else:
        
        return 0.
# specificity (Sp): correclty classified as absent
def Sp(x,y):
    
    if float(sum(x+y==0) + sum(x-y==1))!=0:
        
        return float(sum(x+y==0)) / float(sum(x+y==0) + sum(x-y == 1))  

    else:
        
        return 0.
        
# True skill statistic
def TSS(x,y):
    
    return Sp(x,y)+Sn(x,y)-1
    
# Jaccard index  
def Jaccard(x,y):
    
    x = np.array(x,dtype=float)
    y = np.array(y,dtype=float)
    
    return float(sum(x+y==2))/float(sum(x-y==1)+sum(x-y==-1)+sum(x+y==2))

#  adjusted Jaccard index: penalty on over- and underprediction
def aJaccard(x,y):
    
    x = np.array(x,dtype=float)
    y = np.array(y,dtype=float)
    
    return float(sum(x+y==2)-sum(x-y==1)-sum(x-y==-1))/float(sum(x-y==1)+sum(x-y==-1)+sum(x+y==2))    

# RMSE
def fRMSE(x,y):
    
    return np.sqrt(np.sum((x-y)**2)/float(len(x)))

#
def SSE(x,y):
    
    return np.sum((x-y)**2)

#def wSSE(x,y,w):
#    
#    w_i = np.ones(len(y))
#    w_i[y==1] = w*float(len(y))/float(len(y[y==1]))
#    w_i[y==0] = (1-w)*float(len(y))/float(len(y[y==0]))
#    
#    return np.sum(w_i*(x-y)**2)
    
def nAIC(L,N,K):
    
    vAIC = -2.*np.log(float(L)) + 2.*float(K)
    vAICc = vAIC + float(2*K*(K+1))/float(N-K-1)
    
    return vAIC,vAICc

def nBIC(L,N,K):
    
    #-2.*np.log(float(L))
    BIC = np.log(N)*float(K)
    
    return BIC    

