"""
Created on Tue Feb 26 15:38:35 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import numpy as np
import sys

#def mahalanobis_distance(x,y,normalise = 1.):
#    
#    d  = np.array(x-y)
#    distance = np.sqrt(np.dot(np.dot(np.transpose(d),np.linalg.inv(np.array(cov))),d))
#    
#    if normalise == 1:
#        d = np.zeros(len(x))-np.ones(len(x))
#        reference_distance = np.sqrt(np.dot(np.dot(np.transpose(d),np.linalg.inv(np.array(cov))),d))
#        HSI = distance/reference_distance
#    else:
#        HSI = distance
#        
#    return HSI
    

class EFM():
    """ Environmental filter model used to map n-dimensional hypervolmume and 
    species response
    
    Attributes
    -----------
        'model' (pandas df): coupled environmental-biological input data and 
        model parameters used to compute SI and HSI.
        'nos' (int): number of variables in model
        'flogit' (boolean): flag for logistic increasing and decreasing function
    """
    
    def __init__(self,env_input,model_parameters,logit=False):
        """
        initialise model
        
        Parameters
        ----------
            'env_input' (pandas df): coupled environmental-biological input data
            'model_parameters' (pandas df): model parameters for species response curves        
            'logit' (boolean): flag for logistic increasing and decreasing function
            
        Returns
        -------
            none
            
        """
        self.model = env_input[["sample","ID","X","Y","date","taxon","variable","value"]].merge(model_parameters[["variable","a1","a2","a3","a4","type"]],on="variable",how="inner")
        # number of variables
        number_of_variables = len(self.model["variable"].unique())
        self.nos = number_of_variables
        self.model["value"] = self.model["value"].astype(float)
        self.flogit = logit

    def save_model(self,resname,mode="w"):
        """
        save model to disk
        
        Parameters
        ----------
            'resname' (string): path/dir to which results have to be written to
            'mode' (string): 'w' for write, 'a' for append
            
        Returns
        -------
            none
            
        """        
        if mode=="w":
            self.model.to_csv(resname, mode=mode, header=True)
        else:
            self.model.to_csv(resname, mode=mode, header=False)
    
    def get_model(self):
        """
        return model
        
        Parameters
        ----------
            none
            
        Returns
        -------
            'model' (pandas df): coupled environmental-biological input data and 
            model parameters used to compute SI and HSI.
            
        """          
        return self.model

    def run_models(self):
        """run model
        
        Parameters
        ----------
            none
            
        Returns
        -------
            'model' (pandas df): coupled environmental-biological input data and 
            model parameters used to compute SI and HSI.
            
        """                 
        un_types = self.model["type"].unique().tolist()
        self.model["HSI"]=np.nan

        if self.flogit==False:
            
            for i in un_types:
                
                cond = self.model["type"]==i
    
                if i == "continuous":
                    self.linear(cond)

        else:
            self.logit()
            
        self.model = self.model[~self.model["HSI"].isnull()]
        
        return self.model
            
    def linear(self,condition):
        """
        run linear model (linear decreasing and increasing functions)
        
        Parameters
        ----------
            'condition' (array of bools): records to consider for linear model
            
        Returns
        -------
            none
            
        """                 
        # initiate with HSI == np.nan
        self.model["HSI"]=np.nan
        
        # rule 1: HSI(Value<a1)
        cond = (self.model["value"]<self.model["a1"]) & condition
        self.model.loc[cond,"HSI"] = 0.
        
        # rule 2: HSI(a1<=Value & Value<a2)    
        cond=(self.model["a1"]<=self.model["value"]) & (self.model["value"]<self.model["a2"]) & condition
        self.model.loc[cond,"HSI"] = (self.model["value"][cond]-self.model["a1"][cond])/(self.model["a2"][cond]-self.model["a1"][cond])
        
        # rule 2 exception: HSI(a1=a2=Value)
        cond=(self.model["a1"]==self.model["value"]) & (self.model["value"]==self.model["a2"]) & condition
        self.model.loc[cond,"HSI"] = 1.   
        
        # rule 3: HSI(a2<=value & a3<=value)
        cond=(self.model["a2"]<=self.model["value"]) & (self.model["value"]<=self.model["a3"]) & condition
        self.model.loc[cond,"HSI"] = 1.
        
        # rule 4: HSI(a3<value & a4<=value)
        cond=(self.model["a3"]<self.model["value"]) & (self.model["value"]<=self.model["a4"]) & condition
        self.model.loc[cond,"HSI"]=(self.model["a4"][cond]-self.model["value"][cond])/(self.model["a4"][cond]-self.model["a3"][cond])  
        
        # rule 4 exception: HSI(a3=a4=Value)
        cond=(self.model["value"]==self.model["a3"]) & (self.model["value"]==self.model["a4"]) & condition
        self.model.loc[cond,"HSI"]=1.
        
        # rule 5: HSI(a4<Value)
        cond=(self.model["value"]>self.model["a4"]) & condition
        self.model.loc[cond,"HSI"]=0.

    def logit(self,e=10**-3):
        """
        run logit model (linear decreasing and increasing functions)
        
        Parameters
        ----------
            'e' (float):  minimum defined SI
            
        Returns
        -------
            none
            
        """              
        try:
            self.model["beta1"] = np.log((2-e)/e)/(self.model["a2"]-self.model["a1"]+10**-6)
            self.model["alpha1"] = -self.model["beta1"]*self.model["a2"]
            self.model["beta2"] = np.log((2-e)/e)/(self.model["a3"]-self.model["a4"]+10**-6)
            self.model["alpha2"] = -self.model["beta2"]*self.model["a3"]        
        except:
            print("[PROGRAMMED EXIT] \n\t"+"oops, something went horribly wrong, contact me at sachagobeyn@gmail.com")
            sys.exit("="*19)   
                
        # initiate with HSI == np.nan
        self.model["HSI"]=np.nan
        
        # rule 1: HSI(Value<a1)
        cond = (self.model["value"]<self.model["a1"])
        self.model.loc[cond,"HSI"] = 0.
        
        # rule 2: HSI(a1<=Value & Value<a2)    (logit)
        cond=(self.model["a1"]<=self.model["value"]) & (self.model["value"]<self.model["a2"])
        self.model.loc[cond,"HSI"] = 2/(1+np.exp(-(self.model["beta1"][cond]*self.model["value"][cond]+self.model["alpha1"][cond])))
        
        # rule 3: HSI(a2<=value & a3<=value)
        cond=(self.model["a2"]<=self.model["value"]) & (self.model["value"]<=self.model["a3"])
        self.model.loc[cond,"HSI"] = 1.
        
        # rule 4: HSI(a3<value & a4<=value)
        cond=(self.model["a3"]<self.model["value"]) & (self.model["value"]<=self.model["a4"])
        self.model.loc[cond,"HSI"] = 2/(1+np.exp(-(self.model["beta2"][cond]*self.model["value"][cond]+self.model["alpha2"][cond])))
        
        # rule 5: HSI(a4<Value)
        cond=(self.model["value"]>self.model["a4"])
        self.model.loc[cond,"HSI"] = 0.
        
    def interference(self,mode):
        """interference
    
        Parameters
        ----------
            'mode' (string): type of interference, minimum, mean, squaredproduct or euclidean

        Return
        ------
            aggregated HSI per sample in self.model dataframe
        """
        self.model =  self.model[~self.model["HSI"].isnull()]
        
        if ~self.model.empty:
            
            if mode == "minimum":
                
                return self.model.groupby(["sample","X","Y","date","taxon"]).aggregate({"HSI":np.min}).reset_index()
                
            if mode == "product":
                
                print("[PROGRAMMED EXIT] \n\t"+"implement product formula")
                sys.exit("="*19) 
                
            if mode == "mean":
                
                return self.model.groupby(["sample","X","Y","date","taxon"]).aggregate({"HSI":np.mean}).reset_index()
                
            if mode == "squaredproduct":

                return self.model.groupby(["sample","X","Y","date","taxon"]).aggregate({"HSI":lambda x:(np.nanprod(x))**(1./len(x))}).reset_index()

            if mode == "euclidean":
                    
                return self.model.groupby(["sample","X","Y","date","taxon"]).aggregate({"HSI":lambda x: 1-np.sqrt(np.sum((1.-x)**2))/np.sqrt(len(x))}).reset_index()
                   
        else:
            
            return np.nan
             
#    def threshold(self,threshold):
#        
#        if (~np.isnan(threshold)) & (~np.isnan(self.HSI)):
#        
#            self.present = 0 if self.HSI<threshold else 1
#            
#        else:
#            
#            self.present = np.nan
#            
#        return self.present
#        
#        
#    def diagnose(self):
#        
#        # check if HSI is calculated
#        if ~self.model.empty:
#            
#            return self.model#[self.model["HSI"]!=1.]
#                
#        else:
#            message = "Model is empty"
#            sys.exit(message)
            
#    def project_to_raster(self,raster,X,Y):
#        
#        raster[(X==self.X) & (Y==self.Y)] = self.HSI
#        
#        return raster
def euclidean_distance(x,y,normalise=1):
    """euclidean interference

    Parameters
    ----------
        'x','y' (numpy array): predictions and reference
        'normalize' (int): scale HSI to range 0 - 1
    Return
    ------
        'HSI' (numpy array): aggregated HSI via euclidean distance per sample 
    """    
    
    d  = np.array(x-y)
    distance = np.sqrt(np.sum(d**2))
    #distance = np.sqrt(np.dot(np.dot(np.transpose(d),np.linalg.inv(np.array(cov))),d))
    if normalise == 1:
        d = np.zeros(len(x))-np.ones(len(x))
        #reference_distance = np.sqrt(np.dot(np.dot(np.transpose(d),np.linalg.inv(np.array(cov))),d))
        reference_distance = np.sqrt(np.sum(d**2))
        
        HSI = distance/reference_distance
    else:
        HSI = distance
    return HSI
      
        
