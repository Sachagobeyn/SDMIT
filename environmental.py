"""
Created on Tue Feb 26 15:38:35 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import pandas as pd
import numpy as np
import sys

def environmental_filter_model(X,Y,initial,env_input,species_parameters,model_parameters):

    #initial = initial*mask if mask!=None else initial
    
    " Initialize output "
    output = np.array(initial.copy(),dtype=float)
    
    " Points "
    points = np.unique(env_input["ID"])
    
    for i in range(0,len(points)):
        
        " initiate calculation point "
        pixel = Pixel(env_input["X"][env_input["ID"]==points[i]],env_input["Y"][env_input["ID"]==points[i]])
        
        " contruct models "
        pixel.construct_model(env_input,species_parameters["environmental_parameters"])

        " calculate model according to shape"
        pixel.run_model()
        
        " interference "
        pixel.interference(model_parameters["interference"],threshold=0.5)
        
        " project on raster "
        output = pixel.project_to_raster(output,X,Y)
        
    return output

def mahalanobis_distance(x,y,normalise = 1.):
    
    d  = np.array(x-y)
    distance = np.sqrt(np.dot(np.dot(np.transpose(d),np.linalg.inv(np.array(cov))),d))
    
    if normalise == 1:
        d = np.zeros(len(x))-np.ones(len(x))
        reference_distance = np.sqrt(np.dot(np.dot(np.transpose(d),np.linalg.inv(np.array(cov))),d))
        HSI = distance/reference_distance
    else:
        HSI = distance
        
    return HSI
    
def euclidean_distance(x,y,normalise = 1.0):
        
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
      
class Pixel():
    
    def __init__(self,X=0.,Y=0.):
        
        # we should really fix these lines of code
        try:
            self.X = X.tolist()[0]
            self.Y = Y.tolist()[0]
        except:
            self.X = X
            self.Y = Y
                        
        self.HSI = np.nan
        
    def save_model(self,resname,mode="w"):
        
        if mode=="w":
            self.model.to_csv(resname, mode=mode, header=True)
        else:
            self.model.to_csv(resname, mode=mode, header=False)
    
    def get_model(self):
        
        return self.model
        
    def construct_model(self,env_input,model_parameters):
        
        self.model = env_input[["variable","value"]].merge(model_parameters[["variable","a1","a2","a3","a4","type"]],on="variable",how="right")
        # number of variables
        number_of_variables = len(self.model["variable"].unique())
        self.nos = number_of_variables
        
        return self.model
        
    def run_models(self):
        
        un_types = self.model["type"].unique().tolist()
        self.model["HSI"]=np.nan

        for i in un_types:
            
            cond = self.model["type"]==i

            if i == "acute":
                self.acute(cond)
            if i == "right":
                self.right(cond)
            if i == "line":
                self.line(cond)
        
        return self.model
        
    def line_model(self,condition):
        
        # negative skewed model      
        cond = (self.model["value"]<=self.model["a1"]) & condition
        self.model["HSI"][cond] = 1
        cond = (self.model["value"]>self.model["a1"]) & condition
        self.model["HSI"][cond] = 0
        
        # positive skewed model           
#        cond = (self.model["value"]<self.model["a1"]) & condition
#        self.model["HSI"][cond] = 0
#        cond = (self.model["value"]>=self.model["a1"]) & condition
#        self.model["HSI"][cond] = 1
        
    def right(self,condition):
        
        # 1. negative skewed model      cute.right
        # 1.1. optimal
        cond = (self.model["value"]<=self.model["a1"]) & condition
        self.model["HSI"][cond] = 1
        # 1.2. sub-optimal
        cond = (self.model["a1"]<self.model["value"]) & (self.model["value"]<self.model["a2"]) & condition
        self.model["HSI"][cond] = (self.model["a2"][cond] - self.model["value"][cond])/(self.model["a2"][cond]-self.model["a1"][cond]) 
        # 1.3. non-optimal
        cond = (self.model["value"]>self.model["a2"]) & condition 
        self.model["HSI"][cond] = 0
#        
#        # 2. positive skewed model           
#        cond_model = (model["side"]==1)
#        # 2.1. optimal
#        cond = (model["value"]>=model["a2"]) & cond_model
#        model["HSI"][cond] = 1
#        # 2.2. sub-optimal
#        cond = (model["a1"]<model["value"]) & (model["value"]<model["a2"]) & cond_model
#        model["HSI"][cond] = (model["value"] - model["a1"])/(model["a2"]-model["a1"])     
#        # 2.3. non-optimal
#        cond = (model["value"]<model["a1"]) & cond_model
#        model["HSI"][cond] = 0
#        
    def acute(self,condition):
        
        # initiate with HSI == np.nan
        self.model["HSI"]=np.nan
        
        # rule 1: HSI(Value<a1)
        cond = (self.model["value"]<self.model["a1"]) & condition
        self.model["HSI"][cond] = 0.
        
        # rule 2: HSI(a1<=Value & Value<a2)    
        cond=(self.model["a1"]<=self.model["value"]) & (self.model["value"]<self.model["a2"]) & condition
        self.model["HSI"][cond] = (self.model["value"][cond]-self.model["a1"][cond])/(self.model["a2"][cond]-self.model["a1"][cond])
        
        # rule 2 exception: HSI(a1=a2=Value)
        cond=(self.model["a1"]==self.model["value"]) & (self.model["value"]==self.model["a2"]) & condition
        self.model["HSI"][cond] = 1.   
        
        # rule 3: HSI(a2<=value & a3<=value)
        cond=(self.model["a2"]<=self.model["value"]) & (self.model["value"]<=self.model["a3"]) & condition
        self.model["HSI"][cond] = 1.
        
        # rule 4: HSI(a3<value & a4<=value)
        cond=(self.model["a3"]<self.model["value"]) & (self.model["value"]<=self.model["a4"]) & condition
        self.model["HSI"][cond]=(self.model["a4"][cond]-self.model["value"][cond])/(self.model["a4"][cond]-self.model["a3"][cond])  
        
        # rule 4 exception: HSI(a3=a4=Value)
        cond=(self.model["value"]==self.model["a3"]) & (self.model["value"]==self.model["a4"]) & condition
        self.model["HSI"][cond]=1.
        
        # rule 5: HSI(a4<Value)
        cond=(self.model["value"]>self.model["a4"]) & condition
        self.model["HSI"][cond]=0.

    def interference(self,mode,cov=pd.DataFrame(),threshold=np.nan):
        """interference
    
        Keyword arguments:
        mode -- the name of the file with all model parameters
        Optional arguments:
        cov -- dataframe with covariances (index = variable names, columns= variable names) 
        """
        self.model =  self.model[~self.model["HSI"].isnull()]
        
        if ~self.model.empty:
            
            if mode == "minimum":
                
                self.HSI = float(np.min(self.model["HSI"]))
                
            if mode == "product":
                
                self.HSI = float(np.product(np.array(self.model["HSI"][~self.model["HSI"].isnull()])))
                
            if mode == "mean":
                
                self.HSI = float(np.mean(np.array(self.model["HSI"])))
                
            if mode == "squared product":

                self.HSI = (float(np.prod(np.array(self.model["HSI"]))))**(1./self.nos)
            
            if mode == "mahalanobis":
            
                if cov.empty:
                    message = "Please provide a covariance matrix"
                    sys.exit(message)
                else:
                    # prepare
                    var = self.model["variable"].unique().tolist()                    
                    cov = cov.reset_index()
                    self.model = self.model.merge(cov,on="variable",how="left")
                    self.HSI = mahalanobis_distance(self.model["HSI"],np.ones(len(var)))

            if mode == "euclidean":
                    
                    var = self.model["variable"].unique().tolist() 
                    # if one of the variables, out of range, not present
                    if float(np.min(self.model["HSI"]))==0.:
                        self.HSI = 0.
                    else:
                        self.HSI = euclidean_distance(self.model["HSI"],np.ones(len(var)))
                   
        else:
            
            self.HSI = np.nan
            
        return self.HSI
    
             
    def threshold(self,threshold):
        
        if (~np.isnan(threshold)) & (~np.isnan(self.HSI)):
        
            self.present = 0 if self.HSI<threshold else 1
            
        else:
            
            self.present = np.nan
            
        return self.present
        
        
    def diagnose(self):
        
        # check if HSI is calculated
        if ~self.model.empty:
            
            return self.model#[self.model["HSI"]!=1.]
                
        else:
            message = "Model is empty"
            sys.exit(message)
            
    def project_to_raster(self,raster,X,Y):
        
        raster[(X==self.X) & (Y==self.Y)] = self.HSI
        
        return raster

        