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
      (for example "python script.py -inputdata 'inputdata.csv' -variables 'considered_variables.csv' -taxon 'Baetidae' -resmap 'Results'"
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
      
      fold                     use data either for model development (development) or optimisation (optimisation)
      
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
        model_settings["a1"] = 0
        model_settings["a2"] = 25 
        model_settings["a3"] = 75
        model_settings["a4"] = 100
    
        return model_settings
        
    "Genetic algorithm parameters"
    if filetype=="ga":
        
        ga_settings = {}
        
        ga_settings["number_of_chromosomes"] = 100 
        ga_settings["selection_rate"] = 0.5
        ga_settings["crossover_rate"] = 0.8
        ga_settings["mutation_rate"] = 0.005
        ga_settings["objective_function"] = "AICc"        
        
        "additional"
        ga_settings["maximum_runs"] = 100
        ga_settings["stop_criteria"] = 30
        ga_settings["ncpu"] = -1 #multiprocessor with all available processors
        ga_settings["mode"] = "binary"
        ga_settings["full_output"] = False
        ga_settings["nan_value"] = 100000
    
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
                settings[line[0]] = bool(line[1])
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
    inputdata.to_csv(os.path.join(resmap,"inputdata.csv"))
    #inputdata = pd.read_csv("inputdata.csv")
    
    "Environmental filter parameter estimation (EFPE)"
    from EFPE import EFPE
    model_parameters = EFPE(inputdata,model_settings,taxon,resmap)

    "Optimise model"
    optimisation(inputdata,model_parameters,taxon,variables,ga_settings,resmap)
    #gridcalibration(inputdata,model_parameters,taxon,variables,ga_settings,resmap)
        
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
    model_inputs["taxon"] = taxon
    model_inputs["variables"] = variables
    
    "Initiate boundary conditions for genetic algorithm"
    columns = ["parameter"]
    boundaries = pd.DataFrame(data=np.zeros([len(variables),1]),columns=columns)
    boundaries["parameter"]=variables
    boundaries["cond"]=""
    boundaries["sample"]=1.
    boundaries["down"]=np.nan
    boundaries["up"]=np.nan    

    from GA import GA
    GA("SDM",model_inputs,boundaries,ga_settings,resmap)

def gridcalibration(inputdata,parameters,taxon,variables,ga_settings,resmap):

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
    model_inputs["variables"] = variables
    model_inputs["taxon"] = taxon
    
    "Initiate boundary conditions for genetic algorithm"
    columns = ["parameter"]
    boundaries = pd.DataFrame(data=np.zeros([len(variables),1]),columns=columns)
    boundaries["parameter"]=variables
    boundaries["cond"]=""
    boundaries["sample"]=1.
    boundaries["down"]=np.nan
    boundaries["up"]=np.nan 
 
    from itertools import product
    import datetime
    
    "Create grid"
    grid = list(product(range(2),repeat=len(boundaries)))
    
    "Open time"
    f = open(os.path.join(resmap,"time-grid.txt"),"w")   
    my_dt_ob = datetime.datetime.now()
    time = [my_dt_ob.month, my_dt_ob.day, my_dt_ob.hour, my_dt_ob.minute, my_dt_ob.second]
    f.write("%i,%i,%i,%i,%i\n"%(time[0],time[1],time[2],time[3],time[4]))
    f.close()
    
    "Multiprocessing"
    import multiprocessing
    ncpu = ga_settings["ncpu"]
        
    # get number of processors
    if ncpu==-1:
        ncpu=multiprocessing.cpu_count()
    pool=multiprocessing.Pool(ncpu)
    jobs=[0.]*len(grid)
    
    for i in range(len(grid)):
        jobs[i] = pool.apply_async(eval("gridpoint"),(grid[i],model_inputs,boundaries,ga_settings,resmap,i))  
    pool.close()

    for i in range(len(grid)):
        p = jobs[i].get()
    pool.join()

    f = open(os.path.join(resmap,"time-grid.txt"),"a")    
    my_dt_ob = datetime.datetime.now()
    time = [my_dt_ob.year, my_dt_ob.month, my_dt_ob.day, my_dt_ob.hour, my_dt_ob.minute, my_dt_ob.second]
    f.write("%i,%i,%i,%i,%i\n"%(time[0],time[1],time[2],time[3],time[4]))
    f.close()
    
def gridpoint(point,model_inputs,boundaries,ga_settings,resmap,i) :

    from GA import Chromosoom
    import datetime

    print(str(i)+":"+str(datetime.datetime.now()))

    boundaries.loc[:,'sample'] = point
    chromosoom_i = Chromosoom(boundaries)
    interference(model_inputs,boundaries,chromosoom_i,ga_settings["nan_value"],resmap,False)
            
    return 1
    
def run_filter_model(model_input,taxon,variables,resmap):  
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
    inputdata = inputdata[inputdata["variable"].isin(variables)]  
    inputdata = inputdata[inputdata["fold"]=="optimisation"]
    parameters = model_input["parameters"]
    parameters = parameters[parameters["variable"].isin(variables)]
    
    "Run environmental model"
    from environmental import Pixel
    
    ids = inputdata["ID"].unique().tolist();ind = 0
    
    "Run the model over IDs"
    for j in ids:
        
        # append to output
        cond = (inputdata["ID"]==j)
        test = inputdata[cond]   
        
        model = Pixel(0.,0.)
        model.construct_model(test,parameters)
        model.run_models()
        out = model.get_model()
        # delete nan
        out = out[~out["HSI"].isnull()]
        out["ID"] = j
        out["taxon"] = taxon
        if ind == 0:
            output = deepcopy(out)
        else:
            output = output.append(out)
    
        ind +=1
        
    " Write output "
    output.to_csv(os.path.join(resmap,"SI_values.csv"))
    
    return output    

def SDM(model_input,boundaries,chromosoom,nan_value,res,full_output):    

    """ 
    Function for interference suitability indices
    
    Arguments:
        'model_input' (dictionary): 
            'inputdata' (pandas df): Biological and environmental measurements
                                columns: ["ID","taxon","abundance","variable","value",optional="development"]  
            'parameters' (pandas df): Estimated parameters for habitat preference curves
                                columns: ["taxon","a1","a2","a3","a4","type"]   
        'boundaries' (pandas df): Dataframe with defined boundaries for optimisation problem
        			colums: ["parameter","cond",sample","down","up"]
        			"sample": 0/1 (not) included
        			"parameter": parameter to optimise (in this case filter)
        			"cond": condition applied to sample
        			"down": lower boundary sample
        			"up": upper boundary sample
	'chromosoom' (object):  Chromosoom object (see GA.py)
	'nan_value' (float): value for nan_values
	'res' (str): name of map to write output
	'full_output' (boolean): write full output (yes = 1/no = 0)
    
    Returns:
        'criteria': dictionary with performance indices
 
    """
    
    variables  = mapper(model_input,chromosoom)
    model_output = run_filter_model(model_input,model_input["taxon"],variables,res)
    model_input = model_input["data"]

    "Check whether model has variables"
    "IF model has not variables: set evaluation criterion very high ELSE: run model"
    if len(variables)==0:

        criteria = {}
        criteria = {c:nan_value for c in ["Sn","Sp","TSS","AUC","SSE","AICc","Kappa","CCI"]}

    else:
        
        "Inteference"
        #interference_formula = "lambda x:1-np.sqrt(np.sum((1-x)**2))/np.sqrt(float(len(x)))"
        interference_formula = "lambda x: (np.prod(x))**(1./len(x))"
        model_output = model_output.groupby(["ID","taxon"]).aggregate({"HSI":eval(interference_formula)}).reset_index()
    
        "Link observed abundance to model output"
        model_output = model_output.merge(model_input[["ID","taxon","abundance"]][model_input["fold"]=="optimisation"],how="left",on=["ID","taxon"]).drop_duplicates()
        
        "Evaluate model output on presence/absence"
        model_output["abundance"][model_output["abundance"]!=0] = 1
        
        "Write output"
        if full_output==True:
            model_output.to_csv(os.path.join(res,str(chromosoom.ID)+"-sim.csv"))
        from performance import calculate_performance
        criteria = calculate_performance(model_output,variables,nan_value)
        
    return criteria
    
def mapper(model_input,chromosoom):
    
    "Translate chromosoom to variable inclusion"
    variables = chromosoom.parameters["parameter"][chromosoom.parameters["sample"]==1].tolist()

    return variables
 
def create_dir(res,L):
    
    for i in range(len(L)):
        if not os.path.exists(os.path.join(res,L[i])):
            os.makedirs(os.path.join(res,L[i]))            

if __name__ =="__main__":
    
    "Read parameterfile"
    #arguments = read_parameterfile(sys.argv[1])
    arguments =  read_parameterfile("../parameterfile.txt")
    
    "Overwrite with system parameters"
    arguments = overwrite_arguments(arguments)    

    "Get settings defined in the model and ga setting file"
    model_settings = read_settingsfile(arguments["md_settings"],filetype="md")
    ga_settings = read_settingsfile(arguments["ga_settings"],filetype="ga")
    ga_settings["full_output"] = False
    "Run code"
    run(arguments["inputdata"],arguments["taxon"],arguments["variables"],arguments["resmap"],model_settings,ga_settings)
