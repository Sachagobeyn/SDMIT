# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:34:52 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""
""" 
------------
INTRODUCTION
------------ 
The species distribution model identification tool is a software tool that aims 
to identify (optimise) species distribution models with genetic algorithms. The
implemented SDM is a habitat suitability model, aiming to define the relation 
of species response (here presence/absence) as a function of environmental 
gradients. The genetic algorithm serves as an optimisation method performing:
    
    (1) Input variable selection: Searching for a set of near-optimal input 
                                  parameters which best describe observed pres-
                                  ence/absence patterns (parameters of species
                                  response curves have to be defined a priori).
                                  
    (2) Parameter estimation (PE):  Searching for a set of near-optimal species 
                                    response parameters which best describe 
                                    observed presence/absence patterns.
                                    
    (3) Model identification (MI): Combination of (1) and (2)

The species response curves are defined by trapezoid curves or logit implemen-
tations (logit).

The genetic algorithm encoding has different implemententations for (1), (2) 
    and (3):
    
    (1) Binary: a string of 0 and 1s with a length equal to the  number of 
                    cosidered input variables.
    (2) Continuous: a string of real-values (Haupt and Haupt, 2004) bounded by 
                        predefined boundaries b. Length is equal to the number 
                        of input variables multiplied by 4.
    (3) List of list:   a list of list with real-values bounded by predefined 
                        boundaries b. The first order list is a binary string,
                        whereas a second order continuous string is defined 
                        when a "1" is defined in the first order list. IVS is 
                        encoded in the first order list, PE in the second order

Single objective or multi-objective optimisation is possible, for single object-
ive optimisation True Skill Statistic, Kappa, Sum of Squared Errors, .. can be 
used. For definition of these measures see Mouton et al. (2010). Multi-objecti-
ve optimisation is also feasable, by searching for a Pareto frontier with a non
-dominated sorting genetic algorithme approach (Deb, 2002). Any of the imple-
mented objective can be used.

For implementation of the genetic operators: see
    ADD

For settings of values for hyper parameters, one is advised to follow Gibbs et
    al. (2008) as they where found to be near-optimal.
        
        (1) Determine function evaluations by dividing the computer time 
            available by the average time to compute the fitness function
        (2) Solve equation 1.1 to determine the number of chromosomes:
                
            FE/N log (1-1/N) = -M-log(sqrt(l/12))           equation 1.1
            
            with M = 3 and l = (maximum) length of chromosomes
            for variable length encoding = compute maximum length
            
        (3) Compute mutation rate by dividing 5 by number of chromosomes
        (4) Use elist approach and set crossover rate to 1
            (selection rate = 0.5)

References:
-----------
Deb, K., Pratap, A., Agarwal, S., Meyarivan, T., 2002. A fast and elitist multi-
    objective genetic algorithm: NSGA-II. IEEE Trans. Evol. Comput. 6, 182–197.
Gibbs, M.S., Dandy, G.C., Maier, H.R., 2008. A genetic algorithm calibration 
    method based on convergence due to genetic drift. Inf. Sci. (Ny). 
    178, 2857–2869.
Haupt, R.L., Haupt, S.E., 2004. Algorithms Practical Genetic Algorithms, 2nd 
    ed. John Wiley & Sons, Inc., Hoboken.
Mouton, A.M., De Baets, B., Goethals, P.L.M., 2010. Ecological relevance of 
    performance criteria for species distribution models. Ecol. Modell. 221, 
    1995–2002.
   

------------
INSTRUCTIONS
------------ 

      SDM optimisation is run by running 'script.py' in Python 
      OR in command line 'python script.py -flags'
      
      Implemented for batch simulations
      
      The "parameterfile.txt"-file defines  the input for the code whereas the 
      "settings.txt"-file specify the settings for the model and genetic algorithm
      
      -----------------
      parameterfile.txt
      -----------------
      
      inpudata [STRING]     .csv file with inputdata 
                            (col: [ID,X,Y,taxon,date,abundance,value,variable])

      taxon [STRING]        name of taxon to develop and optimise model
                            (make sure name is in column 'taxon' in inputdata)  
 
      variables [STRING]    .csv file with variables to consider in model
                            (col: [variable,consider,name_sim], name_sim are 
                            the variables checked in inputdata, beware captions,
                            special characters, ..)
                            
      filter_parameters [STRING]  .csv file with parameters of species response 
                                  curves
                            
      resmap [STRING]       name of map to print results
 
      settings [STRING]     name of file containing options for genetic algo-
                            rithm and SDM.
                            
      full_ouput [True OR   True or False
                  False]

      ------------
      inpudata.csv
      ------------
      list with columns (ID,taxon,abundance,variable,value)
      
      ID [INT]              unique ID of point
      
      X [FLOAT]             X coordinate of point, can include nans
      
      Y [FLOAT]             Y coordinate of point, can include nans
              
      date [~]              date, can include nans
          
      taxon [STRING]        name of observed taxon
      
      abundance [INT        presence/absence or abundance of taxon
                 OR FLOAT]
      
      variable  [STRING]    name of measured variable (names must be compatible
                            with name_sim in variables file in parameterfile.txt)
      
      value [FLOAT]         value of the environmental variable
            
      -------------
      variables.csv
      -------------
      list with columns (variable,consider,name_sim)
      
      name_sim [STRING]     name of variable in model (names must be compatible
                            inputdata.csv)
      
      variables [STRING]    synonym or explanation of variable
      
      consider [BINARY]     consider variable in optimisation

      -------------
      filter_parameters.csv
      -------------
      list with columns (taxon,value,type,b1->b4,low,high,a1->a4,variable)
      
      variable [STRING]     name of variable in model (names must be compatible
                            inputdata.csv)
      
      taxon [STRING]        name of observed taxon

      a1 [FLOAT]            lower boundary of range (SI!=0) of species response 
                            curve, bounded by [low,b2]
      
      a2 [FLOAT]            lower boundary of optimal range (SI = 1) of species 
                            response curve, bounded by [b1,b3]

      a3 [FLOAT]            upper boundary of optimal range (SI = 1) of species 
                            response curve, bounded by [b2,b4]

      a4  [FLOAT]           upper boundary of range (SI != 0) of species resp-
                            ponse curve, bounded by [b3,high]
      
      low, high, b0, b1,    boundaries for a (the boundaries b are adjusted                                      
      b2, b3 and b4         dynamicaly when best solution in genetic algorithm
      [FLOAT]               goes out of bound)
      
      --------------
      settings.txt
      --------------
      
      nan_value [FLOAT]             value for nan's in optimisation 
                                    (default: 100000)
      
      number_of_chromosomes [INT]   number of chromosomes to run GA
                                    (default: 100)
    
      selection_rate [FLOAT]        selection rate of population per generation
                                    (default: 0.5 -> 50%)
                                    
      crossover_rate [FLOAT]        crossover rate of population
                                    (default: 1. -> 100%)
      
      mutation_rate [FLOAT]         mutation rate of genes in chromosomes
                                    (default: 0.05 -> 5%)
      
      multi_objective [True         'True' or 'False'
                      OR False]     (default: 'False')
     
      objective_function [STRING]   the objective function to minimize
                                    (default: TSS)
     
      maximum_runs [INT]            maximum number of generations to calculate
                                    (default: 200)
                                         
      ncpu [INT]                    number of cpu's to use on machine 
                                    (-1 = all)
     
      mode [STRING]                 'binary' (IVS), 'continuous' (PE), variable
                                    (MI) (default: MI)

      full_output [True 
                  OR False]         print all outputs of every chromosoom
                                    (default: 'False')
     
    
      ------
      flags (use by running command line python script.py -flags)
      ------
      
      -i                            name of inputdata file
      
      -v                            name of variables file
      
      -t                            name of the taxon
    
      -res                          name of directory to write output
      
      -fp                           parameters of SDM model (species response 
                                     curves)
      -set                          settingsfile
      
      -logit                        True or False
      
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
    keys = ["inputdata","taxon","variables","filter_parameters","resample","resmap","settings","w","full_output","logit"]
    
    "Read lines"
    with open(parameterfile) as f: 
        lines = []
        for line in f:
            lines.append(line.rstrip('\n'))    
    
    "If parameter is defined, overwrite"
    for line in lines:
        
        line = re.split(r'[;,\s]\s*', line)
        
        if line[0] in keys:

            code_arguments[line[0]] = line[1]
            
    return code_arguments

def read_settingsfile(settingfile):    
    """ 
    Get algorithm parameters for model development and optimisation from parameterfile
    
    Arguments:
        'settingsfile' (str): name of file with settings
        'type' (str): either the settingsfile for the model development (md)
                    or model optimisation (ga)
    
    Returns:
       'settings' (dictionary): see function 'default_parameters'
    
    """
    settings = default_settings()
    settings = read_file(settingfile,settings)
    
    return settings
    
def default_settings():
    """ 
    Define default settings for model development and optimisation
    
    Arguments:
        'filetype' (str): either model development (md) or optimisation (ga)
    
    Returns:
       'model_settings' (dictionary): default settings for the model development
       'ga_settings' (dictionary): default settings for the genetic algorithm 
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
  
    #used for run_filter_model, userdefined)   
    settings["interference"] = "squaredproduct"
    settings["logit"] = False
    settings["threshold"] = "max"
    
    return settings

def check_settings(settings):
    
    if (settings["multi_objective"]==True) and (type(settings["objective_function"])==str):
          
          print("[PROGRAMMED EXIT] 'multi_objective = True' but only one objective function is given.")
          sys.exit("="*19)    

    if (settings["multi_objective"]==False) and (type(settings["objective_function"])==list):
          
          print("[PROGRAMMED EXIT] 'multi_objective = False' but two objective functions are given.")
          sys.exit("="*19)   

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

    flags["fp"] = "filter_parameters"
    flags["set"] = "settings"
    flags["logit"] = "logit"
    
    keys = list(flags.keys())
    
    for i in keys:
        
        if "-"+i in sys.argv:
            
            loc = sys.argv.index("-"+i)
            arguments[flags[i]] = eval(sys.argv[loc+1]) if (sys.argv[loc+1]=="True") | (sys.argv[loc+1]=="False") else sys.argv[loc+1] 
        
    return arguments
 
def overwrite_settings(settings):
    """ 
    Overwrite the settings found in the settingsfile file by the system arguments
    
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

def run(inputdata,taxon,variables,filter_parameters,resmap,settings,full_output=False):
    
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
    
    "Load and sample data"
    from data_processing import load_and_preproces_data
    inputdata,filter_parameters,variables = load_and_preproces_data(inputdata,taxon,filter_parameters,variables,resmap,settings["nan_value"])

    "Check input for the optimisation"
    check_input_code(inputdata,taxon,filter_parameters,variables)

    "Optimise model"
    performance,chromosomes = optimisation(inputdata,taxon,filter_parameters,variables,settings,resmap,full_output=full_output)

    "Print performance"
    print_performance(performance,resmap)

    "Print models"
    print_models(inputdata,filter_parameters,taxon,settings["interference"],settings,chromosomes,resmap)

def check_input_code(inputdata,taxon,filter_parameters,variables):
    
    cond = False
    
    error = []
    
    if taxon not in inputdata["taxon"].tolist():

        error.append("Taxon in parameterfile is not found inputdata file")
        cond = True

    if taxon not in filter_parameters["taxon"].tolist():

        error.append("Taxon in parameterfile is not found model parameter file")
        cond = True
    
    if cond==False:
        
        if np.sum(inputdata["variable"].isin(variables))==0:

            error.append("None of the variables in variable file not found in inputdata file")
            cond = True
        
        if np.sum(filter_parameters["variable"].isin(variables))==0:

            error.append("None of the variables in variable file not found in filter parameter file")
            cond = True
   
    if cond == True:
    
        print("[PROGRAMMED EXIT] \n\t"+"\n \t".join(error))
        sys.exit("="*19) 
  
def optimisation(inputdata,taxon,filter_parameters,variables,settings,resmap,full_output=False):

    """ 
    Function to run model and optimise the developed model.
    
    Arguments:
        'inputdata' (pandas df): Biological and environmental measurements
                            columns: ["sample","ID","taxon","abundance","variable","value",optional="development"]  
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
    filter_parameters = filter_parameters.sort_values(["variable","value"],ascending=True)
    model_inputs["parameters"] = filter_parameters
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
    boundaries = filter_parameters

    "Check input for model and raise warning if model will be empty"
    check_input_code(inputdata,taxon,filter_parameters,variables)
    
    from GA import GA    
    performance,solution = GA("translate_chromosome2model",model_inputs,boundaries,settings,resmap,full_output=full_output)
    
    return performance,solution
    
def translate_chromosome2model(model_inputs,boundaries,chromosoom,mode,nan_value,resmap,full_output,final_run=False):

    "Step 1: Extract variables which are included in model"
    variables = chromosoom.parameters["variable"][~(chromosoom.parameters["sample"]==0)].tolist()
    
    "Step 2: Extract initial guess for parameters"
    parameters = model_inputs["parameters"][model_inputs["parameters"]["variable"].isin(variables)]
 
    "Step 3: deepcopy inputs for model to redefine model inputs with chromosomes model parameters"
    "(if no deepcopy than model parameters  of other chromosomes are linked)"
    inputs = deepcopy(model_inputs)
    
    "Step 3: Check whether model is not empty"
    if len(variables)!=0:
        
        "Step 3a: Transform chromosome encoding to model parameters"
        inputs["parameters"],inputs["K"] = translate_chromosome(deepcopy(parameters),chromosoom,mode)
        
        "Step 3b: Run the model"
        performance = run_filter_model(inputs,resmap,chromosoom.ID,full_output=full_output)
        
    else:
        performance = {c:nan_value if (c=="AIC") or (c=="SSE") or (c=="BIC") else -nan_value for c in ["AIC","N","K","SSE","AUC","Kappa","CCI","Sp","Sn","TSS","threshold","BIC"]}

    if performance["N"]!=len(model_inputs["data"]["ID"].unique()):
        
        performance = {c:nan_value if (c=="AIC") or (c=="SSE") or (c=="BIC") else -nan_value for c in ["AIC","N","K","SSE","AUC","Kappa","CCI","Sp","Sn","TSS","threshold","BIC"]}

    "Step 4: Print the model if required"
    if full_output==True:
        parameters.to_csv(os.path.join(resmap,str(int(chromosoom.ID))+"-parameters.csv"))
    
    return performance,parameters

def translate_chromosome(parameters,chromosoom,mode):
    
    K=0
                
    if mode=="binary":

        un_var = chromosoom.parameters["variable"][chromosoom.parameters["sample"]==1].tolist()

        parameters =  parameters[parameters["variable"].isin(un_var)]
        
        factor = 1
        K = np.sum(parameters["type"]=="discrete") + np.sum(parameters["type"]=="continuous")*factor
        
    
    if (mode=="variable") or (mode=="continuous"):

        "Overwrite a"
        for i in range(1,5):
            parameters["a"+str(i)] = np.nan
     
        
        un_var = parameters["variable"].unique().tolist()
        chromosoompar = chromosoom.parameters
        
        for i in un_var:
            
            cond_par = parameters["variable"]==i
            cond_chr = chromosoompar["variable"]==i
            string = chromosoompar["sample"][cond_chr].values[0]
            
            "discrete/acute"
            if parameters["type"][cond_par].iloc[0]=="categorical":
                string =  string.returnString()
                parameters.loc[cond_par,["a"+str(i) for i in range(1,len(string)+1,1)]] = string
                K += len(string)
            else:
                parameters.loc[cond_par,["a1","a2","a3","a4"]] = string.returnString()
                K += 4
                             
    return parameters,K
    
def initiate_continues(model_inputs,parameters,type_filter="acute"):

    abiotic_filters = parameters["variable"][parameters["type"]==type_filter].tolist()
    boundaries = pd.DataFrame(data=np.zeros([len(abiotic_filters)*2,1]),columns=["parameter"])
    boundaries["parameter"] = [i+"_"+["a2","a3"][j] for i in abiotic_filters for j in range(2)]
    boundaries["sample"]=[parameters[i.split("_")[1]][parameters["variable"]==i.split("_")[0]].values[0] for i in boundaries["parameter"].unique()]
    boundaries["down_cond"]=[('boundaries["sample"][boundaries["parameter"]=="'+str(i.replace("a3","a2")) +'"].values[0]') if ("a3" in i) else (str(parameters["a1"][parameters["variable"]==i.split("_")[0]].values[0])) for i in boundaries["parameter"].unique()] 
    boundaries["up_cond"]=[('boundaries["smaxample"][boundaries["parameter"]=="'+str(i.replace("a2","a3")) +'"].values[0]') if ("a2" in i) else (str(parameters["a4"][parameters["variable"]==i.split("_")[0]].values[0])) for i in boundaries["parameter"].unique()]
    boundaries["range"] = [parameters["a3"][parameters["variable"]==i.split("_")[0]].values[0]-parameters["a2"][parameters["variable"]==i.split("_")[0]].values[0] for i in boundaries["parameter"].unique()] 
    boundaries["C"] = 1
    
    return boundaries
    
def initiate_discrete(model_inputs,parameters,type_filter="discrete"):
    
    boundaries = deepcopy(parameters[["variable","value","a1"]][parameters["type"]==type_filter])
    boundaries["parameter"] = boundaries["variable"]+"_"+boundaries["value"]
    boundaries["sample"] = boundaries["a1"]
    boundaries["down_cond"] = "0"
    boundaries["up_cond"] = "1"
    boundaries["range"] = 1
    boundaries["C"] = 0
    
    return boundaries[["parameter","sample","down_cond","up_cond","range","C"]]
    
def run_filter_model(model_input,resmap,ID,full_output=False):
    
    """ 
    Function to run environmental filter model
    
    Arguments:
        'model_input' (dictionary): 
            'inputdata' (pandas df): Biological and environmental measurements
                                columns: ["sample","ID","taxon","abundance","variable","value",optional="development"]  
            'parameters' (pandas df): Estimated parameters for habitat preference curves
                                columns: ["taxon","a1","a2","a3","a4","type"]   
        'taxon' (str): name of taxon                        
        'variables' (list): self-explanatory
        'resmap' (str): name of map to write output
        
    
    Returns:
        'output' (pandas df): output of environmental filter model
                            columns: ["sample","ID","taxon","variable","value","a1","a2","a3","a4","type","HSI"]
 

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
    output["prediction"] = output["HSI"]*output["RSI"]
    
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

    else:
        
        threshold = float(threshold)
        
    perf = calculate_performance(output,K,evaluate=True,threshold=threshold)  
    
    return perf     
   
def print_optimised_model(parameters,opt_parameters,resmap):
    
    variables = []
    opt_parameters = opt_parameters.reset_index()
    
    for i in range(len(opt_parameters)):
        
        variable,par = opt_parameters["parameter"].iloc[i].split("_")
        parameters[par][parameters["variable"]==variable] = opt_parameters["sample"].iloc[i]
        variables.append(variable)
    
    return parameters[parameters["variable"].isin(variables)]
        
def create_dir(res,L):
    
    for i in range(len(L)):
        if not os.path.exists(os.path.join(res,L[i])):
            os.makedirs(os.path.join(res,L[i]))            

def print_performance(performance,resmap):
    
    f = open(os.path.join(resmap,"performance.csv"),"w")
    f.write("Criterion,Value\n")
    for i in performance.keys():
        f.write(i+","+str(performance[i])+'\n')
    f.close()
    
def print_models(inputdata,filter_parameters,taxon,interference,settings,chromosomes,resmap):
    
    model_input = {}
    model_input["data"]= inputdata
    model_input["interference"]=interference
    model_input["settings"] = settings
    model_input["logit"] = settings["logit"]
    model_input["threshold"] = settings["threshold"]
    for i in chromosomes:
        variable = i.parameters["variable"][~(i.parameters["sample"]==0)].tolist()
        par_i,K = translate_chromosome(deepcopy(filter_parameters[filter_parameters["variable"].isin(variable)]),i,settings["mode"])
#        cond =  ~par_i["grid"].isnull()
#        par_i["grid"][cond] = par_i["grid"][cond].apply(lambda x:x.returnGrid())
        par_i.to_csv(os.path.join(resmap,"optimisation_summary","parameters_"+str(i.ID)+".csv"))
        if len(par_i)>0:
            model_input["parameters"] = par_i
            model_input["K"] = K
            #model_input["Kmax"] = 
            run_filter_model(model_input,os.path.join(resmap,"model_runs"),i.ID,full_output=True)

    variables = chromosomes[0].parameters["variable"][~(chromosomes[0].parameters["sample"]==0)].tolist()
    opt_parameters,_ = translate_chromosome(deepcopy(filter_parameters.loc[filter_parameters["variable"].isin(variables)]),chromosomes[0],mode=settings["mode"])
    opt_parameters.to_csv(os.path.join(resmap,"optimal_parameters_"+taxon+".csv"))
    
if __name__ =="__main__":
 
    print("="*19)
    print("SDMIT (version 2)")
    print("="*19)    
    print(" CC BY 4.0 Creative Commons, sacha gobeyn, \n"+
          " sacha.gobeyn at ugent.be OR sachagobeyn at gmail.com")
    print(" https://github.com/Sachagobeyn/SDMIT \n")
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
    run(arguments["inputdata"],arguments["taxon"],arguments["variables"],arguments["filter_parameters"],arguments["resmap"],settings,full_output=eval(arguments["full_output"]))

