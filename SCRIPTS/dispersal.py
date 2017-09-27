"""
Created on Tue Feb 26 15:38:35 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import pandas as pd
import numpy as np
import sys

class DFM():
    
    def __init__(self,env_input,model_parameters):
                    
        # stack dataframes
        for i in range(4):
            
            parameters_i = pd.melt(model_parameters,id_vars=['taxon','variable','b'+str(i)], var_name=['none'],value_vars='a'+str(i+1));
            parameters_i["a1"] = parameters_i["value"]
            parameters_i["value"] = parameters_i["b"+str(i)]
            parameters_i = parameters_i[["taxon","variable","value","a1"]]
            if i==0:
                parameters = parameters_i  
            else:
                parameters = parameters.append(parameters_i)
            
        model_parameters = parameters[~parameters["value"].isnull()]
        
        self.model = env_input[["sample","ID","X","Y","date","taxon","variable","value"]].merge(model_parameters[["variable","value","a1"]],on=["variable","value"],how="inner")
        # number of variables
        number_of_variables = len(self.model["variable"].unique())
        self.nos = number_of_variables
        
    def save_model(self,resname,mode="w"):
        
        if mode=="w":
            self.model.to_csv(resname, mode=mode, header=True)
        else:
            self.model.to_csv(resname, mode=mode, header=False)
    
    def get_model(self):
        
        return self.model

    def run_models(self):
        
        self.model["RSI"]=np.array(self.model["a1"])
        self.model = self.model[~self.model["RSI"].isnull()]
        
    def interference(self,mode):
        """interference
    
        Keyword arguments:
        mode -- the name of the file with all model parameters
        Optional arguments:
        cov -- dataframe with covariances (index = variable names, columns= variable names) 
        """
        self.model =  self.model[~self.model["RSI"].isnull()]
        
        if ~self.model.empty:
            
            if mode == "minimum":
                
                print("implement")
                
            if mode == "product":
                
                print("implement")
                
            if mode == "mean":
                
                return self.model.groupby(["sample","X","Y","date","taxon"]).aggregate({"RSI":np.mean}).reset_index()
                
            if mode == "squared product":

                return self.model.groupby(["sample","X","Y","date","taxon"]).aggregate({"RSI":lambda x:(np.prod(x))**(1./len(x))}).reset_index()

            if mode == "euclidean":
                    
                return self.model.groupby(["sample","X","Y","date","taxon"]).aggregate({"RSI":lambda x: 1-np.sqrt(np.sum((1.-x)**2))/np.sqrt(len(x))}).reset_index()
                   
        else:
            
            return np.nan
        
        
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

        
