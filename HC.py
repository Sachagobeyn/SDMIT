# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 10:11:07 2016

@author: sacha
"""

import numpy as np
import operator
import copy
import random
from script import *
import shutil
import os
import pandas as pd
pd.set_option('chained_assignment',None)
from runtime import *
from itertools import product

def HC(model,model_inputs,boundaries,HC_options,resmap,full_output=False):
    
    print("="*21)
    print("Start HC optimisation")
    print("="*21)
    
    # create directory
    if os.path.exists(os.path.join(resmap,"HC_results")):
        shutil.rmtree(os.path.join(resmap,"HC_results"))    
    if os.path.exists(os.path.join(resmap,"HC_model_parameters")):
        shutil.rmtree(os.path.join(resmap,"HC_model_parameters"))    
    create_dir(resmap,["HC_results","HC_model_parameters"])
#    
#    f = open(os.path.join(res,"model_runs","cpu.txt"),"w") 
#    f.close()

    # extract model options
    neighbourhood = HC_options["neighbourhood"]
    objective_function = "SSE"
    init_of = HC_options["init_of"]
    full_output = HC_options["full_output"]
    
    if "ncpu"  in HC_options:
	ncpu = int(HC_options["ncpu"])
    else:
	ncpu = -1
 
    vc = Neighbour(boundaries)
    vc.setEvaluationcriteria(init_of)
    
    run = 0; cond = 0
    
    while cond<5:
        
        # Step 1: Initialise vn
        neighbours = initialise_vn(vc,neighbourhood)
        
        # Step 2: Run neighbours
        neighbours,dt = multithread(neighbours,model,model_inputs,objective_function,os.path.join(resmap,"HC_model_parameters"),ncpu=ncpu)#,boundaries,objective_function,nan_value,os.path.join(res,"model_runs"),ncpu,full_output=full_output)
        
        # Step 3: Get best neighbour
        best = evaluate_neighbours(neighbours)

        # Step 4: Print neighbours
        if full_output==True: 
            print_neighbours(neighbours,os.path.join(resmap,"HC_results"),run)        
        
        # Step 4: Compare with vc        
        if (best.evaluation_criteria<vc.evaluation_criteria):
            vc = deepcopy(best)
        else:
            cond += 1
            neighbourhood += 0.05
        run += 1

        print("Iteration ("+str(run)+"): best = "+str(vc.evaluation_criteria) +"\t (runtime: "+str(dt)+" seconds)")
    
    print_neighbours(neighbours,os.path.join(resmap,"HC_results"),"results_ga")       
    print("="*18)
    print("Stop HC optimisation")
    print("="*18)
    
    return best           

class Neighbour():
    
    def __init__(self,parameters):
        
        self.parameters = parameters
        self.performance = np.nan
        self.ID = np.nan
        self.evaluation_criteria = np.nan
        
    def setPerformance(self,performance):
        
        self.performance = performance
    
    def setEvaluationcriteria(self,evaluation_criteria):
        
        self.evaluation_criteria = evaluation_criteria
        
    def setID(self,n):
        
        self.ID = n
        
def initialise_vn(vc,neighbourhood):
    
    #it_n = list(product([0.,0.],repeat=len(vc.parameters)))
    it_n = list(product([-neighbourhood,0.,neighbourhood],repeat=len(vc.parameters)))
    n = 0
    
    neighbours = []
    
    for i in range(len(it_n)):
        
        parameters_i = deepcopy(vc.parameters)
        parameters_i["sample"] =  vc.parameters["sample"]+vc.parameters["range"]*np.array(it_n[i])
        parameters_i = check_boundaries(parameters_i)
        neighbour_i = Neighbour(parameters_i)
        neighbour_i.ID = n;n+=1
        neighbours.append(neighbour_i)
        
    return neighbours
    
def check_boundaries(boundaries):
    
    boundaries["down"] = np.nan
    boundaries["up"] = np.nan
    
    for i in range(len(boundaries)):
        
        boundaries["down"].iloc[i] = eval(boundaries["down_cond"].iloc[i])
        boundaries["up"].iloc[i] = eval(boundaries["up_cond"].iloc[i])
        
    boundaries["sample"][boundaries["sample"]<boundaries["down"]]=boundaries["down"]
    boundaries["sample"][boundaries["sample"]>boundaries["up"]]=boundaries["up"]
    
    return boundaries
    
def multithread(neighbours,model,model_input,objective_function,resmap,ncpu=-1):
    
    # import package
    import multiprocessing

    # start timer
    nos = np.sum([1. if np.isnan(neighbours[i].evaluation_criteria) else 0. for i in range(len(neighbours))])
    t = Runtime(nos,10)
    
    # get number of processors
    if ncpu==-1:
        ncpu=multiprocessing.cpu_count()

    # pool processors
    pool=multiprocessing.Pool(ncpu)

    
    # make list of jobs
    jobs=[0.]*len(neighbours)
    for i in range(len(neighbours)):
        if np.isnan(neighbours[i].evaluation_criteria):
            jobs[i] = pool.apply_async(eval(model),(model_input,neighbours[i],resmap))  
    pool.close()
    
    for i in range(len(neighbours)):
        if np.isnan(neighbours[i].evaluation_criteria):  
            performance = jobs[i].get()
            neighbours[i].setPerformance(performance)
            neighbours[i].setEvaluationcriteria(performance[objective_function])
    
    pool.join()
    dt = t.close(print_value=False)
    
    return neighbours,dt     

def evaluate_neighbours(neighbours):
    
    # sort chromosomes on performance
    neighbours.sort(key=operator.attrgetter('evaluation_criteria'),reverse=False)
    
    return neighbours[0]

def print_neighbours(neighbours,res,run):

    # performance
    perf_keys = neighbours[0].performance.keys();
    # columns    
    columns = list(neighbours[0].parameters["parameter"]) + ["ID"] + perf_keys
    
    # data
    data = np.zeros([len(neighbours),len(columns)])
    
    for i in range(len(neighbours)):
        
        data[i,:-(1+len(perf_keys))] = np.array(neighbours[i].parameters["sample"]).ravel()
        data[i,-(1+len(perf_keys))]  = neighbours[i].ID
        data[i,-len(perf_keys):]  = np.array([neighbours[i].performance[perf_keys[j]] for j in range(len(perf_keys))])
    
    dataframe = pd.DataFrame(data=data,columns=columns)
    #pd.DataFrame(data=data,columns=columns).to_csv(os.path.join(res,str(run)+".csv"))
    dataframe.to_csv(os.path.join(res,str(run)+".csv"))
    