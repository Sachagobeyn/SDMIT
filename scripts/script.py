# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:34:52 2015
Description: see https://github.com/Sachagobeyn/SDMIT/ or https://github.com/Sachagobeyn/SDMIT/releases/
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""

import pandas as pd
import numpy as np
import os
from copy import deepcopy
import re 
import sys
#pd.set_option('chained_assignment',None)
import warnings
warnings.filterwarnings("error")

cwd = os.getcwd()

def read_parameterfile(parameterfile):

    """ get code arguments from the parameterfile
    
    Parameters
    ----------        
    'parameterfile' (str): name of the parameterfile
    
    Returns
    -------
    'code_arguments' (dictionary):
        
            'inputdata' (str): name of inputdata 
            'variables' (str): name of variables
            'taxon' (str): name of the taxon
            'resmap' (str): name of directory to write output
            'model_parameters' (str): name of model_parameters
            'full_output' (bool): write out all model runs
            'logit' (bool): use logistic increasing and decreasing function to
                            describe species response
    """
    code_arguments = {}
    keys = ["inputdata","taxon","variables","model_parameters","resmap","settings","full_output","logit"]
    
    "Read lines for parameterfile"
    with open(parameterfile) as f: 
        lines = []
        for line in f:
            lines.append(line.rstrip('\n'))    
    
    "If parameter is defined in the parameterfile, overwrite default"
    for line in lines:
        
        line = re.split(r'[;,\s]\s*', line)
        
        if line[0] in keys:

            code_arguments[line[0]] = line[1]
    
    if 'scripts' in cwd:
        
        for i in ["inputdata","variables","model_parameters","resmap","settings"]:
            
            code_arguments[i] = os.path.join("..",code_arguments[i])

    return code_arguments

def read_settingsfile(settingfile):    
    """ get optimisation settings from parameterfile
    
    Parameters
    ----------  
        'settingsfile' (str): name of file with settings
        
    Returns
    ----------  
       'settings' (dictionary): for explanation settings see function
       'default_settings'
    
    """
    settings = default_settings()
    settings = read_file(settingfile,settings)
    
    return settings
    
def default_settings():
    """ define settings for optimisation
    
    Parameters
    ----------  
        none
    
    Returns
    -------
       'settings' (dictionary): default settings for running optimisation
       
    """
     
    settings = {}
    
    "Options for SDM"
    settings = default_settings_structure(settings)

    " Settings for evoluationary algorithm"
    settings["number_of_chromosomes"] = 100 
    settings["selection_rate"] = 0.5
    settings["crossover_rate"] = 1.
    settings["mutation_rate"] = 0.05
    settings["maximum_runs"] = 100
    settings["stop_criteria"] = 30
    
    # adaptive hyper parameter control was tested but omitted from final code 
    settings["adaptive"] = False
    settings["k1"] = 1.
    settings["k2"] = 0.5
    settings["k3"] = 1.
    settings["k4"] = 0.5
    
    " Allow duplicate chromosomes in evolving solution space"
    settings["duplicates"] = False
    
    " Objective function "
    settings["multi_objective"] = False
    settings["objective_function"] = "TSS"

    " Subject to optimisation, input variable selection, parameter estimation or both"
    settings["mode"] = "binary" #mode = binary for IVS, variable for IVS + PE

    " Settings for multiprocessing and nan_values in evolutionary algorithm"
    settings["ncpu"] = -1 #multiprocessor with all available processors
    #settings["precision"] = 3 #bit precison of problem (discretization)
    settings["nan_value"] = 100000
    
    return settings
   
def default_settings_structure(settings):
    """ define default settings for optimisation
    
    Parameters
    ----------  
        'settings' (dictionary): default settings for the genetic algorithm 
    
    Returns
    -------
       'settings' (dictionary): default settings for the genetic algorithm 
                                updated with default SDM settings
       
    """
     
    
    #used for run_filter_model, userdefined)   
    settings["interference"] = "squaredproduct"
    settings["logit"] = False
    settings["threshold"] = "max"
    
    return settings

def check_settings(settings):
 
    """ check settings for optimisation
    
    Parameters
    ----------  
        'settings' (dictionary): default settings for optimisation
    
    Returns
    -------
       none
       
    """
    
    if (settings["multi_objective"]==True) and (type(settings["objective_function"])==str):
          
          print("[PROGRAMMED EXIT] 'multi_objective = True' but only one objective function is given.")
          sys.exit("="*19)    

    if (settings["multi_objective"]==False) and (type(settings["objective_function"])==list):
          
          print("[PROGRAMMED EXIT] 'multi_objective = False' but two objective functions are given.")
          sys.exit("="*19)   

def read_file(setting_file,settings):
    """ read settingsfile line by line
    
    Parameters
    ----------  
        'setting_file' (string): name settingsfile
    
    Returns
    -------
         'settings' (dictionary): settings for optimisation   
    """    
    
    "Read lines"
    with open(setting_file) as f: 
        lines = []
        for line in f:
            lines.append(line.rstrip('\n'))
    
    "If parameter is defined, overwrite"
    for line in lines:
        
        line = re.split(r'[;,\s]\s*', line)
        if line[0] in settings:
            "Check type"
            "String/float"
            if (len(line)>2) and (line[2]!=""):
                settings[line[0]] =[line[i] for i in range(1,len(line))]
            else:
                if isinstance(settings[line[0]],str):
                    settings[line[0]] = str(line[1])
                elif isinstance(settings[line[0]],bool):
                    settings[line[0]] = eval(line[1])
                else:
                    settings[line[0]] = float(line[1])
                

    return settings

def overwrite_arguments(arguments):
    """ overwrite parameterfile arguments with system arguments (flags)
    
    Parameters
    ----------  
        'arguments' (dictionary): parameter file arguments
    
    Returns
    -------
        'arguments' (dictionary): updated parameter file arguments
        NOTE: flags are identifiers recognised by the system by a dash (-)
    """    
        
    flags = {}
    flags["i"] = "inputdata"
    flags["t"] = "taxon"
    flags["res"] = "resmap"
    flags["v"] = "variables"

    flags["fp"] = "model_parameters"
    flags["set"] = "settings"
    flags["logit"] = "logit"
    
    keys = list(flags.keys())
    
    for i in keys:
        
        if "-"+i in sys.argv:
            
            loc = sys.argv.index("-"+i)
            arguments[flags[i]] = eval(sys.argv[loc+1]) if (sys.argv[loc+1]=="True") | (sys.argv[loc+1]=="False") else sys.argv[loc+1] 
        
    return arguments
 
def overwrite_settings(settings):
    """ overwrite settings with system arguments (flags)
    
    Parameters
    ----------  
        'settings' (dictionary): settings for optimisation 
    
    Returns
    -------
        'settings' (dictionary): updated settings for optimisation 
        NOTE: flags are identifiers recognised by the system by a dash (-)
    """  
    
    flags = {}
    flags["adaptive"] = "adaptive"
    flags["duplicates"] = "duplicates"
    flags["local_optimisation"] = "local_optimisation"
    flags["pc"] = "crossover_rate"
    flags["pm"] = "mutation_rate"
    flags["k1"] = "k1"
    flags["k2"] = "k2"
    flags["thr"] = "threshold"

    keys = list(flags.keys())
    
    for i in keys:
        
        if "-"+i in sys.argv:
            
            loc = sys.argv.index("-"+i)
            settings[flags[i]] = eval(sys.argv[loc+1]) if (sys.argv[loc+1]=="True") | (sys.argv[loc+1]=="False") else sys.argv[loc+1] 
            
    return settings

def run(inputdata,taxon,variables,model_parameters,resmap,settings,full_output=False):
    """ overwrite settings with system arguments (flags)
    
    Parameters
    ----------   
        'inputdata' (pandas df): input data
        'taxon' (str): name of taxon
        'variables' (pandas df): considered variables 
        'model_parameters' (str): model parameters for species response curves        
        'resmap' (str): name of directory to write output
        'settings' (dictionary): settings for optimisation         
        'full_output' (bool): write out all model runs
        
    Returns
    -------
        none        
    """
    
    "Create output map"
    create_dir("",[resmap])
    
    "Load and sample data"
    from data_processing import load_and_preproces_data
    inputdata,model_parameters,variables = load_and_preproces_data(inputdata,taxon,model_parameters,variables,resmap,settings["nan_value"])

    "Check input for the optimisation"
    check_input_code(inputdata,taxon,model_parameters,variables)

    "Optimise model"
    performance,chromosomes = optimisation(inputdata,taxon,model_parameters,variables,settings,resmap,full_output=full_output)

    "Print performance"
    print_performance(performance,resmap)

    "Print models"
    print_models(inputdata,model_parameters,taxon,settings,chromosomes,resmap)

def check_input_code(inputdata,taxon,model_parameters,variables):
    """ check for input arguments to optimisation function
    
    Parameters
    ----------   
        'inputdata' (pandas df): input data
        'taxon' (str): name of taxon
        'model_parameters' (str): model parameters for species response curves        
        'variables' (pandas df): considered variables 
        
    Returns
    -------
        none        
    """    
    cond = False
    
    error = []
    
    if taxon not in inputdata["taxon"].tolist():

        error.append("Taxon in parameterfile is not found inputdata file")
        cond = True

    if taxon not in model_parameters["taxon"].tolist():

        error.append("Taxon in parameterfile is not found model parameter file")
        cond = True
    
    if cond==False:
        
        if np.sum(inputdata["variable"].isin(variables))==0:

            error.append("None of the variables in variable file not found in inputdata file")
            cond = True
        
        if np.sum(model_parameters["variable"].isin(variables))==0:

            error.append("None of the variables in variable file not found in filter parameter file")
            cond = True
   
    if cond == True:
    
        print("[PROGRAMMED EXIT] \n\t"+"\n \t".join(error))
        sys.exit("="*19) 
  
def optimisation(inputdata,taxon,model_parameters,variables,settings,resmap,full_output=False):
    """ main function to initialise optimisation 
    
    Parameters
    ----------   
        'inputdata' (pandas df): input data
        'taxon' (str): name of taxon
        'model_parameters' (str): model parameters for species response curves        
        'variables' (pandas df): considered variables 
         'resmap' (str): name of directory to write output
        'settings' (dictionary): settings for optimisation         
        'full_output' (bool): write out all model runs
        
    Returns
    -------
        'performance' (dictionary): values for evaluation measures, number of 
        data samples, threshold
        'solution' (list): GA.chromosome objects of last iteration, each 
        containing fitness and model parameters
    """    
  
    "Prepare inputs to run model and run the filter model"
    model_inputs = {}
    model_inputs["data"]= inputdata
    model_parameters = model_parameters.sort_values(["variable","value"],ascending=True)
    model_inputs["parameters"] = model_parameters
    model_inputs["interference"]=settings["interference"]
    model_inputs["settings"] = settings
    model_inputs["logit"] = settings["logit"]
    model_inputs["threshold"] = settings["threshold"]
    
    "Initiate boundary conditions for genetic algorithm"
    columns = ["parameter"]
    boundaries = pd.DataFrame(data=np.zeros([len(variables),1]),columns=columns)
    boundaries["parameter"]=variables
    boundaries["cond"]=""
    boundaries["sample"]=1.
    
    "Define boundary for model"
    boundaries = model_parameters

    "Check input for model and raise warning if model will be empty"
    check_input_code(inputdata,taxon,model_parameters,variables)
    
    from GA import GA    
    performance,solution = GA("compute_fitness",model_inputs,boundaries,settings,resmap,full_output=full_output)
    
    return performance,solution
    
def compute_fitness(model_inputs,boundaries,chromosome,mode,nan_value,resmap,full_output,final_run=False):
    """ function which computes fitness
    
    Parameters
    ----------   
        'model_inputs' (dictionary):  contains inputdata, interference flag, 
        logit flag, model parameters for species response curves, threshold flag,
        number of input variables K
        'boundaries' (pandas df): dataframe containing model boundary conditions 
        for specific taxon (format: model_parameters df)
        'chromosome' (GA.chromosome object): object containing fitness and candidate 
        model parameters
        'mode' (string): 'variable' (embedded) or 'binary' (wrapper) feature
        selection
        'nan_value' (float): nan value for computation
        'resmap' (string): path to results map
        'full_output' (boolean)
        'final_run': flag for indicating if iteration is last cycle
        
    Returns
    -------
        'performance' (dictionary): values for evaluation measures, number of 
        data samples, threshold
        'parameters' (pandas df): parameters for species response curve
    """    

    "Step 1: Extract variables which are included in model"
    variables = chromosome.parameters["variable"][~(chromosome.parameters["sample"]==0)].tolist()
    
    "Step 2: Extract initial guess for parameters"
    parameters = model_inputs["parameters"][model_inputs["parameters"]["variable"].isin(variables)]
 
    "Step 3: deepcopy inputs for model to redefine model inputs with chromosomes model parameters"
    "(if no deepcopy than model parameters  of other chromosomes are linked)"
    inputs = deepcopy(model_inputs)
    
    "Step 3: Check whether model is not empty"
    if len(variables)!=0:
        
        "Step 3a: Transform chromosome encoding to model parameters"
        inputs["parameters"],inputs["K"] = translate_chromosome(deepcopy(parameters),chromosome,mode)
        
        "Step 3b: Run the model"
        performance = run_model(inputs,resmap,chromosome.ID,full_output=full_output)
        
    else:
        performance = {c:nan_value if (c=="AIC") or (c=="SSE") or (c=="BIC") else -nan_value for c in ["AIC","N","K","SSE","AUC","Kappa","CCI","Sp","Sn","TSS","threshold","BIC"]}

    if performance["N"]!=len(model_inputs["data"]["ID"].unique()):
        
        performance = {c:nan_value if (c=="AIC") or (c=="SSE") or (c=="BIC") else -nan_value for c in ["AIC","N","K","SSE","AUC","Kappa","CCI","Sp","Sn","TSS","threshold","BIC"]}

    "Step 4: Print the model if required"
    if full_output==True:
        parameters.to_csv(os.path.join(resmap,str(int(chromosome.ID))+"-parameters.csv"))
    
    return performance,parameters

def translate_chromosome(parameters,chromosome,mode):
    """ function which translates genotype encoded in chromosome to model with species response curve parameters
    
    Parameters
    ----------   
        'parameters' (pandas df): initial parameters for species response curve
        'chromosome' (GA.chromosome object): object containing fitness and model 
        parameters
        'mode' (string): 'variable' (embedded) or 'binary' (wrapper) feature
        selection
        
    Returns
    -------
        'parameters' (pandas df): candidate parameters for species response curve
        extracted from chromosome
        'K': number of input variables K
    """      
    K=0
                
    if mode=="binary":

        un_var = chromosome.parameters["variable"][chromosome.parameters["sample"]==1].tolist()

        parameters =  parameters[parameters["variable"].isin(un_var)]
        
        factor = 1
        K = np.sum(parameters["type"]=="discrete") + np.sum(parameters["type"]=="continuous")*factor
        
    
    if (mode=="variable") or (mode=="continuous"):

        "Overwrite a"
        for i in range(1,5):
            parameters["a"+str(i)] = np.nan
     
        
        un_var = parameters["variable"].unique().tolist()
        chromosoompar = chromosome.parameters
        
        for i in un_var:
            
            cond_par = parameters["variable"]==i
            cond_chr = chromosoompar["variable"]==i
            string = chromosoompar["sample"][cond_chr].values[0]
            
            "categorical/continuous"
            if parameters["type"][cond_par].iloc[0]=="categorical":
                string =  string.returnString()
                parameters.loc[cond_par,["a"+str(i) for i in range(1,len(string)+1,1)]] = string
                K += len(string)
            else:
                parameters.loc[cond_par,["a1","a2","a3","a4"]] = string.returnString()
                K += 4
                             
    return parameters,K
        
def run_model(model_input,resmap,ID,full_output=False):
    """ run model 
    
    Parameters
    ----------   
        'model_inputs' (dictionary):  contains inputdata, interference flag, 
        logit flag, candidate model parameters for species response curves, 
        threshold flag, number of input variables K
        'resmap' (string): path to results map
        'full_output' (boolean)
        'ID' (int): unique tag of chromosome
        
    Returns
    -------
        'performance' (dictionary): values for evaluation measures, number of 
        data samples, threshold
    """ 
    
    "Extract inputdata and habitat preference curve parameters"
    inputdata = model_input["data"]
    interference = model_input["interference"]
    logit = model_input["logit"]
    parameters=model_input["parameters"]
    threshold=model_input["threshold"]
    K =  model_input["K"]
    
    "Run dispersal filter model"
    from dispersal import DFM
    if np.sum(parameters["type"]=="categorical")>0:
        model = DFM(inputdata,parameters[parameters["type"]=="discrete"])
        model.run_models()
        #model.save_model(os.path.join(resmap,str(neighbour)+"_DF.csv"))
        output_DF = model.interference(interference)
    else:
        output_DF = []
        
    "Run environmental model"
    from environmental import EFM
    if np.sum(parameters["type"]=="continuous")>0:
        model = EFM(inputdata,parameters[parameters["type"]=="continuous"],logit=logit)
        model.run_models()
        #model.save_model(os.path.join(resmap,str(neighbour)+"_EF.csv"))
        output_EF = model.interference(interference)
    else: 
        output_EF = []
        
    "Interfere"   
    if len(output_DF)==0:
        output = output_EF
        output["RSI"] = 1.
    if len(output_EF)==0:
        output = output_DF
        output["HSI"] = 1.
    if (len(output_EF)!=0) & (len(output_DF)!=0):
        output = output_DF.merge(output_EF,on=["X","Y","date","sample","taxon"],how="inner")

    "Link observed abundance to model output"
    output = output.merge(inputdata[["X","Y","date","sample","taxon","abundance"]],how="left",on=["sample","taxon"]).drop_duplicates()
    

    "Get prediction presence/absence"
    output.loc[:,"prediction"] = output["HSI"]*output["RSI"]
    
    "Evaluate model output on presence/absence"
    output.loc[:,"observation"] = 0
    output.loc[output["abundance"]!=0,"observation"] = 1
    
    "Print"    
    if full_output==True:

        output.to_csv(os.path.join(resmap,str(int(ID))+"-model_run.csv"))
     
    "Calculate criteria"
    from performance import calculate_performance
    
    if threshold == "max":
        
        threshold = np.nan
        
    elif threshold == "prev":
        
        threshold = 1-np.sum(output["observation"])/len(output)

    elif threshold == "prob":
    
        threshold = "prob"
        
    else:
        
        threshold = float(threshold)
        
    performance = calculate_performance(output,K,evaluate=True,threshold=threshold)  
    
    return performance     
           
def create_dir(resmap,L):
    """ create directory for output to which results are written to
    
    Parameters
    ----------   
        'resmap' (str): name/path of main output directory
        
    Returns
    -------
        'L' (list): list of names which have to be written under res directory
    """ 
    for i in range(len(L)):
        if not os.path.exists(os.path.join(resmap,L[i])):
            os.makedirs(os.path.join(resmap,L[i]))            

def print_performance(performance,resmap):
    """ print performance 'best' model to disk
    
    Parameters
    ----------   
        'performance' (dictionary): values for evaluation measures, number of 
        data samples, threshold of 'best' found model
        'resmap' (str): name/path of main output directory
        
    Returns
    -------
        none
    """  
    f = open(os.path.join(resmap,"performance.csv"),"w")
    f.write("Criterion,Value\n")
    for i in performance.keys():
        f.write(i+","+str(performance[i])+'\n')
    f.close()
    
def print_models(inputdata,model_parameters,taxon,settings,chromosomes,resmap):
    """ print all models of population to disk
    
    Parameters
    ----------   
        'inputdata' (pandas df): input data
        'taxon' (str): name of taxon
        'model_parameters' (str): model parameters for species response curves        
        'variables' (pandas df): considered variables 
        'settings' (dictionary): settings for optimisation
        'chromosomes' (list):   GA.chromosome objects, each 
        containing fitness and model parameters      
         'resmap' (str): name of directory to write output
        
    Returns
    -------
        none
    """      
    model_input = {}
    model_input["data"]= inputdata
    model_input["interference"]=settings['interference']
    model_input["settings"] = settings
    model_input["logit"] = settings["logit"]
    model_input["threshold"] = settings["threshold"]
    for i in chromosomes:
        variable = i.parameters["variable"][~(i.parameters["sample"]==0)].tolist()
        par_i,K = translate_chromosome(deepcopy(model_parameters[model_parameters["variable"].isin(variable)]),i,settings["mode"])
#        cond =  ~par_i["grid"].isnull()
#        par_i["grid"][cond] = par_i["grid"][cond].apply(lambda x:x.returnGrid())
        par_i.to_csv(os.path.join(resmap,"optimisation_summary","parameters_"+str(i.ID)+".csv"))
        if len(par_i)>0:
            model_input["parameters"] = par_i
            model_input["K"] = K
            #model_input["Kmax"] = 
            run_model(model_input,os.path.join(resmap,"model_runs"),i.ID,full_output=True)

    variables = chromosomes[0].parameters["variable"][~(chromosomes[0].parameters["sample"]==0)].tolist()
    opt_parameters,_ = translate_chromosome(deepcopy(model_parameters.loc[model_parameters["variable"].isin(variables)]),chromosomes[0],mode=settings["mode"])
    opt_parameters.to_csv(os.path.join(resmap,"optimal_parameters_"+taxon+".csv"))
    
if __name__ =="__main__":
 
    print("="*19)
    print("SDMIT (version 2)")
    print("="*19)    
    print(" CC BY 4.0 Creative Commons, sacha gobeyn, \n"+
          " sacha.gobeyn at ugent.be OR sachagobeyn at gmail.com")
    print(" https://github.com/Sachagobeyn/SDMIT \n https://github.com/Sachagobeyn/SDMIT/releases/")
    print(" Reading parameterfile and code settings")
       
    "Read parameterfile"
    if len(sys.argv)>1:
        arguments = read_parameterfile(sys.argv[-1])
    else:
        arguments =  read_parameterfile("../parameterfile.txt")
    
    "Overwrite arguments with system parameters"
    arguments = overwrite_arguments(arguments)    

    "Get settings defined in the setting file"
    settings = read_settingsfile(arguments["settings"])
 
    "Overwrite settings with system parameters"
    settings = overwrite_settings(settings) 
    
    "Check settings"
    check_settings(settings)
        
    "Run code"
    run(arguments["inputdata"],arguments["taxon"],arguments["variables"],arguments["model_parameters"],arguments["resmap"],settings,full_output=eval(arguments["full_output"]))

