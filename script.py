# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:34:52 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

""" 
------------
INSTRUCTIONS
------------ 

      Model is run by running 'script.py' in Python OR in command line 'python script.py -flags'
      Implemented for batch simulations
      The "parameterfile.txt"-file is read to the model and optimisation
      
      -----------------
      parameterfile.txt
      -----------------
      
      inpudata                  .csv file with inputdata

      taxon                     name of taxon to develop and optimise model
 
      variables                 .csv file with variables to consider in model

      resmap                    name of map to print results
 
      ga_settings               name of file containing options for genetic optimisation algorithm

      md_settings               name of file containing options for model development
      
      ------------
      inpudata.csv
      ------------
      list with columns (ID,taxon,abundance,variable,value,fold)
      
      ID                       ID of point
      
      taxon                    name of observed taxon
      
      abundance                abundance of taxon
      
      variable                 name of measured variable (names must be compatible with variables file in parameterfile.txt)
      
      value                    value of the variable
      
      fold                     use data either for model development (developmen) or optimisation (optimisation)
      
      -------------
      variables.csv
      -------------
      list with columns (variable,consider)
      
      variable                 name of measured variable (names must be compatible with variables file in parameterfile.txt)
      
      consider                 consider variable in model dev. and opt.
      
      -----------------
      modelsettings.txt
      -----------------
         
      a1,a2,a3,a4               percentiles of enviromental data for parameters of unimodal non-uniform habitat preference curves

      --------------
      gasettings.txt             (for detailed information, see Haupt and Haupt, 2004)
      --------------
      
      number_of_chromosomes     number of chromosomes to run GA
    
      selection_rate            selection rate of population per generation
     
      crossover_rate            crossover rate of population
     
      mutation_rate             mutation rate of genes in chromosomes
     
      objective_function        the objective function to minimize
     
      maximum_runs              maximum number of generations to calculate
     
      stop_criteria             number of generations for which no improvement is observed in the best value of the objective function.
     
      ncpu                      number of cpu's to use on machine (-1 = all)
     
      mode                      binary (a continues varient is implemented but not properly tested)
    
      full_output               print all outputs of every chromosoom tested
     
      nan_value                 nan-value for objective and evalutation criteria for infeasible solutions (has to be a number, int or float). Typically set way above the values of the -to minize!- objective function
    
      ------
      flags
      ------
      
      -inputdata                 name of inputdata file
      
      -variables                 name of variables file
      
      -taxon                     name of the taxon
    
      -resmap                    name of directory to write output
      
"""

import pandas as pd
pd.set_option('chained_assignment',None)
import numpy as np
import os
from copy import deepcopy
import re 
import sys

def read_parameterfile(parameterfile):

    """ 
    Get code arguments from the parameterfile
    
    Arguments:
        'parameterfile' (str): name of the parameterfile
    
    Returns:
        'code_arguments' (dictionary):
        
            'inputdata' (str): name of inputdata 
            NOTE: structure inputdata should be conform Section XX tutorial
            'variables' (str): name of variables
            NOTE: structure variables should be conform Section XX tutorial
            'taxon' (str): name of the taxon
            'resmap' (str): name of directory to write output
    
    """
    code_arguments = {}
    keys = ["inputdata","taxon","variables","resmap","ga_settings","md_settings"]
    
    "Read lines"
    lines = [line.rstrip('\n') for line in open(parameterfile)]
    
    "If parameter is defined, overwrite"
    for line in lines:
        
        line = re.split(r'[;,\s]\s*', line)
        
        if line[0] in keys:

            code_arguments[line[0]] = line[1]
            
    return code_arguments

def read_settingsfile(settingfile,filetype):    
    """ 
    Get algorithm parameters for model development and optimisation from parameterfile
    
    Arguments:
        'settingsfile' (str): name of file with settings
        'type' (str): either the settingsfile for the model development (md)
                    or model optimisation (ga)
    
    Returns:
       'settings' (dictionary): see function 'default_parameters'
    
    """
    settings = default_settings(filetype)
    settings = read_file(settingfile,settings)
    
    return settings
    
def default_settings(filetype):
    """ 
    Define default settings for model development and optimisation
    
    Arguments:
        'filetype' (str): either model development (md) or optimisation (ga)
    
    Returns:
       'model_settings' (dictionary): default settings for the model development
       'ga_settings' (dictionary): default settings for the genetic algorithm 
    """
    
    "Model development parameters (see XXX)"
    if filetype=="md":
        
        model_settings = {}
        model_settings["a1"] = 5
        model_settings["a2"] = 25 
        model_settings["a3"] = 75
        model_settings["a4"] = 95
    
        return model_settings
        
    "Genetic algorithm parameters"
    if filetype=="ga":
        
        ga_settings = {}
        
        ga_settings["number_of_chromosomes"] = 100 
        ga_settings["selection_rate"] = 0.5
        ga_settings["crossover_rate"] = 0.8
        ga_settings["mutation_rate"] = 0.005
        ga_settings["objective_function"] = "AIC"        
        
        "additional"
        ga_settings["maximum_runs"] = 100
        ga_settings["stop_criteria"] = 30
        ga_settings["ncpu"] = -1 #multiprocessor with all available processors
        ga_settings["mode"] = "binary"
        ga_settings["full_output"] = False
        ga_settings["nan_value"] = 100000
        
        "hill climbing"
        ga_settings["neighbourhood"] = 0.05 
    
        return ga_settings
    
def read_file(setting_file,settings):
    """ 
    Read settingsfile
    
    Arguments:
        'setting_file' (str): name of the setting file
        'filetype' (str): either model development (md) or optimisation (ga)
    Returns:
        'settings' (dictionary)
    """
    
    "Read lines"
    lines = [line.rstrip('\n') for line in open(setting_file)]
    
    "If parameter is defined, overwrite"
    for line in lines:
        
        line = re.split(r'[;,\s]\s*', line)
                   
        if line[0] in settings:
            "Check type"
            "String/float"
            if isinstance(settings[line[0]],str):
                settings[line[0]] = str(line[1])
            elif isinstance(settings[line[0]],bool):
                settings[line[0]] = eval(line[1])
            else:
                settings[line[0]] = float(line[1])

    return settings

def overwrite_arguments(arguments):
    """ 
    Overwrite the arguments found in the parameter file by the system arguments
    
    Arguments:
        'arguments' (dictionary):
        
            'inputdata' (str): name of inputdata file
            'variables' (str): name of variables file
            'taxon' (str): name of the taxon
            'resmap' (str): name of directory to write output
            
    Returns:
        'code_arguments' (dictionary):
        
            'inputdata' (str): name of inputdata 
            NOTE: structure inputdata should be conform Section XX tutorial
            'variables' (str): name of variables
            NOTE: structure variables should be conform Section XX tutorial
            'taxon' (str): name of the taxon
            'resmap' (str): name of directory to write output
    
    """
    
    flags = {}
    flags["i"] = "inputdata"
    flags["t"] = "taxon"
    flags["res"] = "resmap"
    flags["v"] = "variables"
    keys = flags.keys()
    
    for i in keys:
        
        if "-"+i in sys.argv:
            
            loc = sys.argv.index("-"+i)
            arguments[flags[i]] = sys.argv[loc+1] 
            
    return arguments
  
def run(inputdata,taxon,variables,resmap,model_settings,ga_settings):
    """ 
    Run script for model development and optimisation
    
    Arguments:
        'inputdata' (str): name of inputdatafile 
        'taxon' (str): name of taxon
        'variables' (str): name of variablesfile
        'resmap' (str): name of map to write output
        'model_setting' (dictionary): settings for model development
        'ga_setting' (dictionary): settings for optimisation (with genetic algorithm)
    
    Returns:

    """
    
    "Create output map"
    create_dir("",[resmap])
    
    "Load and sample data (see XXX)"
    from data_processing import load_and_sample_data
    inputdata, variables = load_and_sample_data(inputdata,taxon,variables,resmap)
    
    "Environmental filter parameter estimation (EFPE)"
    from EFPE import EFPE
    model_parameters = EFPE(inputdata,model_settings,taxon,resmap)

    "Optimise model"
    optimisation(inputdata,model_parameters,taxon,variables,ga_settings,resmap)
    
def optimisation(inputdata,parameters,taxon,variables,ga_settings,resmap):

    """ 
    Function to run model and optimise the developed model.
    
    Arguments:
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["ID","taxon","abundance","variable","value",optional="development"]  
        'parameters' (pandas df): Estimated parameters for habitat preference curves
                            columns: ["taxon","a1","a2","a3","a4","type"]        
        'taxon' (str): name of taxon
        'variables' (pandas df): Lists of considered variables
                            columns: ["variable","consider"]
        'ga_setting' (dictionary): settings for optimisation (with genetic algorithm)                            
        'resmap' (str): name of map to write output
    
    Returns:

    """
  
    "Prepare inputs to run model and run the filter model"
    model_inputs = {}
    model_inputs["data"]= inputdata
    model_inputs["parameters"] = parameters
    model_inputs["output"] = run_filter_model(model_inputs,os.path.join(resmap,"SI_values.csv"))
    
    "Initiate boundary conditions for genetic algorithm"
    columns = ["parameter"]
    boundaries = pd.DataFrame(data=np.zeros([len(variables),1]),columns=columns)
    boundaries["parameter"]=variables
    boundaries["cond"]=""
    boundaries["sample"]=1.
    boundaries["down"]=np.nan
    boundaries["up"]=np.nan    

    from GA import GA    
    best,_ = GA("interference",model_inputs,boundaries,ga_settings,resmap)
    
    "Prepare inputs to run model and run the filter model"
    parameters = parameters[parameters["variable"].isin(best.parameters["parameter"][best.parameters["sample"]==1].tolist())]
    
    model_inputs = {}
    model_inputs["data"]= inputdata
    model_inputs["parameters"] = parameters[parameters["variable"].isin(best.parameters["parameter"][best.parameters["sample"]==1].tolist())]
    
    "Initiate boundary conditions for hill climbing"
    columns = ["parameter"]
    boundaries = pd.DataFrame(data=np.zeros([len(model_inputs["parameters"]["variable"].tolist() )*2,1]),columns=columns)
    boundaries["parameter"] = [i+"_"+["a2","a3"][j] for i in model_inputs["parameters"]["variable"].tolist() for j in range(2)]
    boundaries["sample"]=[model_inputs["parameters"][i.split("_")[1]][model_inputs["parameters"]["variable"]==i.split("_")[0]].values[0] for i in boundaries["parameter"].unique()]
    boundaries["down_cond"]=[('boundaries["sample"][boundaries["parameter"]=="'+str(i.replace("a3","a2")) +'"].values[0]') if ("a3" in i) else (str(model_inputs["parameters"]["a1"][model_inputs["parameters"]["variable"]==i.split("_")[0]].values[0])) for i in boundaries["parameter"].unique()] 
    boundaries["up_cond"]=[('boundaries["sample"][boundaries["parameter"]=="'+str(i.replace("a2","a3")) +'"].values[0]') if ("a2" in i) else (str(model_inputs["parameters"]["a4"][model_inputs["parameters"]["variable"]==i.split("_")[0]].values[0])) for i in boundaries["parameter"].unique()]
    boundaries["range"] = [model_inputs["parameters"]["a4"][model_inputs["parameters"]["variable"]==i.split("_")[0]].values[0]-model_inputs["parameters"]["a1"][model_inputs["parameters"]["variable"]==i.split("_")[0]].values[0] for i in boundaries["parameter"].unique()] 
    
    "Optimise parameters of HPC"
    from HC import HC
    ga_settings["init_of"] = best.evaluation_criteria
    
    best = HC("evaluate_species_parameters",model_inputs,boundaries,ga_settings,resmap)

def run_filter_model(model_input,res,interference=""):  
    """ 
    Function to run environmental filter model
    
    Arguments:
        'model_input' (dictionary): 
            'inputdata' (pandas df): Biological and environmental measurements
                                columns: ["ID","taxon","abundance","variable","value",optional="development"]  
            'parameters' (pandas df): Estimated parameters for habitat preference curves
                                columns: ["taxon","a1","a2","a3","a4","type"]   
        'taxon' (str): name of taxon                        
        'variables' (list): self-explanatory
        'resmap' (str): name of map to write output
        
    
    Returns:
        'output' (pandas df): output of environmental filter model
                            columns: ["ID","taxon","variable","value","a1","a2","a3","a4","type","HSI"]
 

    """
    
    "Extract inputdata and habitat preference curve parameters"
    inputdata = model_input["data"]
    parameters = model_input["parameters"]
    
    "Run environmental model"
    from environmental import EFM
    
    "Run the model over IDs"
    model = EFM(inputdata,parameters)
    model.run_models()
    
    if interference!="":

        output = model.interference(interference)
        "Link observed abundance to model output"
        output = output.merge(inputdata[["ID","taxon","abundance"]],how="left",on=["ID","taxon"]).drop_duplicates()
        "Evaluate model output on presence/absence"
        output["abundance"][output["abundance"]!=0] = 1
        "Calculate criteria"
        output.to_csv(res)
        from performance import calculate_performance
        return calculate_performance(output,model_input["K"],False)      
        
    else:
        
        model.get_model().to_csv(res)
        return model.get_model()

def interference(model_input,boundaries,chromosoom,nan_value,res,full_output):
    
    "Get model input and output"
    model_output = model_input["output"]
    model_input = model_input["data"]
    
    variables = chromosoom.parameters["parameter"][chromosoom.parameters["sample"]==1].tolist()
 
    "Check whether model has variables"
    "IF model has not variables: set evaluation criterion very high ELSE: run model"
    if len(variables)==0:

        criteria = {}
        criteria = {c:nan_value for c in ["AUC","N","K","SSE","AIC","Kappa","CCI","K","Sp","Sn","TSS"]}

    else:
        
        "Consider filters in the candidate model"
        model_output =  model_output[model_output["variable"].isin(variables)]
    
        "Inteference"
        interference_formula = "lambda x:1-np.sqrt(np.sum((1-x)**2))/np.sqrt(float(len(x)))"
        #interference_formula = "lambda x:np.prod(x)**(1./len(x))"
        interference_formula = "np.mean"
        model_output = model_output.groupby(["ID","taxon"]).aggregate({"HSI":eval(interference_formula)}).reset_index()
    
        "Link observed abundance to model output"
        model_output = model_output.merge(model_input[["ID","taxon","abundance"]],how="left",on=["ID","taxon"]).drop_duplicates()
        
        "Evaluate model output on presence/absence"
        model_output["abundance"][model_output["abundance"]!=0] = 1
        
        "Write output"
        if full_output==True:
            model_output.to_csv(os.path.join(res,str(chromosoom.ID)+"-sim.csv"))
        from performance import calculate_performance
        criteria = calculate_performance(model_output,2.*len(variables),evaluate=True)
        
    return criteria

def evaluate_species_parameters(model_input,neighbour,resmap):
    
    parameters = model_input["parameters"]
    inputdata = model_input["data"]
    
    for i in neighbour.parameters["parameter"].tolist():
        parameters[i.split("_")[1]][parameters["variable"]==i.split("_")[0]] = neighbour.parameters["sample"][neighbour.parameters["parameter"]==i].values[0]

    # save parameters 
    model_input = {}
    model_input["data"] = inputdata
    model_input["parameters"] = parameters
    model_input["K"] = len(neighbour.parameters)
    
    return run_filter_model(model_input,os.path.join(resmap,str(neighbour.ID)+".csv"),interference="mean")
    
def create_dir(res,L):
    
    for i in range(len(L)):
        if not os.path.exists(os.path.join(res,L[i])):
            os.makedirs(os.path.join(res,L[i]))            

if __name__ =="__main__":
    
    "Read parameterfile"
    #arguments = read_parameterfile(sys.argv[1])
    arguments =  read_parameterfile("parameterfile.txt")
    
    "Overwrite with system parameters"
    arguments = overwrite_arguments(arguments)    

    "Get settings defined in the model and ga setting file"
    model_settings = read_settingsfile(arguments["md_settings"],filetype="md")
    ga_settings = read_settingsfile(arguments["ga_settings"],filetype="ga")
    
    "Run code"
    run(arguments["inputdata"],arguments["taxon"],arguments["variables"],arguments["resmap"],model_settings,ga_settings)
