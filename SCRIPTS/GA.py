# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 14:30:49 2015
Description:  
@author: sacha gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com)
"""
#def GA(boundary_file,GA_parameters):
#    
#    # load boundary file
#    boundary = pd.read_csv(boundary_file)    
#    
#    # Get GA parameters
#    n_chromosoons = GA["n_chromosoons"]
#        
#    # Step 1: initialize the model parameters
#    chromosoons = [0.]*n_chromosoons
#        

#def sample_parameters(boundaries): 
#    
#    # initialize parameters
#    parameters = boundary.copy()
#    
#    
#    parameters["sample"] =  (boundary["UP"]-boundary["DOWN"])*np.random.rand()+boundary["DOWN"]
#            
#        
#def check_boundaries(sample_parameters,boundary):
#
#    sample_parameters =     
#    
#def read_parameterfile
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
import datetime

def GA(model,model_inputs,boundaries,GA_options,res,full_output=False):
    
    # create directory
    if os.path.exists(os.path.join(res,"calibration_results")):
        shutil.rmtree(os.path.join(res,"calibration_results"))    
    if os.path.exists(os.path.join(res,"model_runs")):
        shutil.rmtree(os.path.join(res,"model_runs"))    
    create_dir(res,["model_runs","calibration_results"])
    
    f = open(os.path.join(res,"model_runs","cpu.txt"),"w") 
    f.close()
    
    # write time in file per generation
    f = open(os.path.join(res,"calibration_results","meta_results.txt"),"w")
    f.write("generation,best,month,day,hour,minute,second\n")
    f.close()

    # extract model options
    n_chromosomes  = int(GA_options["number_of_chromosomes"]);
    selection_rate = GA_options["selection_rate"]
    mutation_rate  = GA_options["mutation_rate"]
    crossover_rate = GA_options["crossover_rate"]
    criteria = GA_options["stop_criteria"]
    objective_function = GA_options["objective_function"]
    mode = GA_options["mode"]
    nan_value = GA_options["nan_value"]
    full_output = GA_options["full_output"]

    if "ncpu"  in GA_options:
        ncpu = int(GA_options["ncpu"])
    else:
        ncpu = -1
 
    # initialize chromosomes and best
    n = 0
    chromosomes,n = initialize_chromosomes(boundaries[["cond","parameter","down","up"]],n_chromosomes,mode,n)
    best = np.nan
    
    # while 
    run = 0 
    best_criteria = []
    cond  = True                
    #t = Runtime(GA_options["runs"],10) 
    
    while cond:
        
        # Stepo 1a: remove duplicates + introduce random chromosomes (keep diversity)
        chromosomes,n = proces_duplicates(chromosomes,'parameters["sample"].tolist()',boundaries,n_chromosomes,n)
        
        # Step 1b: calculate performanc
        #t = time.time()
        chromosomes = multithread(chromosomes,model,model_inputs,boundaries,objective_function,nan_value,os.path.join(res,"model_runs"),ncpu,full_output=full_output)
        #print("One iteration ("+str(n_chromosomes)+"): ")
        #print(datetime.timedelta(seconds=(time.time()-t)))       

        # Step 2a: evaluate chomosomes
        chromosomes,best = evaluate_chromosomes(chromosomes,best)
        
        # Step 2b: print chromosomes
        if run+1==criteria:
            print_chromosooms(chromosomes,os.path.join(res,"calibration_results"),run)   
        
        # Step 3: natural selection (use elitism)
        chromosomes,n_keep = selection(chromosomes,selection_rate)
        
        # Step 4: pairing
        parents = pairing(chromosomes,n_chromosomes,n_keep)
        
        # Step 5: mating
        chromosomes,n = mating(chromosomes,parents,crossover_rate,mode,n_chromosomes,n)
         
        # Step 6: mutations 
        chromosomes,n = mutation(chromosomes,mutation_rate,mode,n)
        
        # Step 7: save and iterate
        ## save best performance 
        best_criteria.append(best.evaluation_criteria)
        
        ## print best to screen
        print("Iteration (%i): best = %.2f"%(run,best.evaluation_criteria))
                
        ## iterate and evaluate while condition
        run+=1
        cond = False if run >= GA_options["maximum_runs"] else True
        print(best_criteria)        
        
        if (run>criteria):

            cond = False if np.sum(np.array(best_criteria[-(int(criteria)):])==best_criteria[-1])==criteria else True
            
        f = open(os.path.join(res,"calibration_results","meta_results.txt"),"a")
        my_dt_ob = datetime.datetime.now()
        time = [my_dt_ob.month, my_dt_ob.day, my_dt_ob.hour, my_dt_ob.minute, my_dt_ob.second]
        f.write("%i,%.2f%i,%i,%i,%i,%i\n"%(run,best.evaluation_criteria,time[0],time[1],time[2],time[3],time[4]))
        f.close()
        #t.iteration(run)
        
    #t.close()
  
    
    
    return best,chromosomes
    
class Chromosoom():
    
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
        
def initialize_chromosomes(boundaries,n_chromosomes,mode,n):
    
    chromosomes = [0.]*n_chromosomes
    
    for i in range(n_chromosomes):
    
        # continues GA
        if mode == "continuous":
            parameters = sample_parameters_continues(boundaries)
        # binary GA (either 0 or 1)
        if mode == "binary":
            parameters = sample_parameters_binary(boundaries)

        chromosomes[i] = Chromosoom(parameters)
        chromosomes[i].setID(n);n+=1

    return chromosomes,n  
      
def sample_parameters_continues(boundaries):
    
    # implement test case
    parameters = boundaries.copy(deep=True)
    parameters["sample"] = np.random.rand(len(parameters))
    
    # check boundaries
    for i in range(len(boundaries)):
        
        if parameters["cond"].loc[i] !="":
            
            cond = parameters["cond"].iloc[i].split(" ")
            
            if cond[0]=="<":

                parameters["sample"].iloc[i]=np.random.uniform(0,parameters["sample"][parameters["parameter"]==cond[1]])    
                
            else:
                
                parameters["sample"].iloc[i]=np.random.uniform(parameters["sample"][parameters["parameter"]==cond[1]],1)
        
    return boundaries[["parameter","sample","down","up","cond"]]

def sample_parameters_binary(boundaries):
    
    parameters = boundaries.copy(deep=True)
    parameters["sample"] = np.random.randint(2.,size=len(parameters))
    
    return parameters

def check_boundary_conditions(parameters):
    
    # fix out of boundar
    parameters["sample"][parameters["sample"]<0] = np.random.rand()
    parameters["sample"][parameters["sample"]>1] = np.random.rand()
    
    # check if conditions apply
    for i in range(len(parameters)):
        
        if parameters["cond"].loc[i] !="":
            
            cond = parameters["cond"].loc[i].split(" ")
            
            if cond[0]=="<":

                if parameters["sample"].loc[i] > parameters["sample"][parameters["parameter"]==cond[1]]:
                    
                    return True
                
            else:
                
                if parameters["sample"].loc[i] < parameters["sample"][parameters["parameter"]==cond[1]]:
        
                    return True
    
        else:
            
            return False

def proces_duplicates(chromosomes,attributes,boundaries,n_chromosomes,n):
    
    chromosomes = remove_duplicates(chromosomes,attributes)
    
    for i in range(n_chromosomes-len(chromosomes)):
        
        new,n = initiate_chromosoom(boundaries,n)
        chromosomes.append(new)
        
    return chromosomes,n
    
def remove_duplicates(objects,attribute):
    
    seen = list()
    unique = []
    for obj in objects:
        if eval('obj.'+str(attribute)) not in seen:
            unique.append(obj)
            seen.append(eval('obj.'+str(attribute)))
    return unique

def initiate_chromosoom(boundaries,n):
    
    # generate n_chromosomes-n_keep random samples
    parameters = sample_parameters_binary(boundaries)
    new = Chromosoom(parameters)
    new.setID(n);n+=1
        
    return new,n

def evaluate_chromosomes(chromosomes,best):
    
    # sort chromosomes on performance
    chromosomes.sort(key=operator.attrgetter('evaluation_criteria'),reverse=False)
    
    # elitism
    if type(best)!=float:
        
        # if there is a chromosome with a better peformanca: replace
        if chromosomes[0].evaluation_criteria > best.evaluation_criteria:
            chromosomes.append(best)
            # sort chromosomes on performance
            chromosomes.sort(key=operator.attrgetter('evaluation_criteria'),reverse=False)
            #print([i.evaluation_criteria for i in chromosomes])
        else:
           # keep best solution
           best = copy.deepcopy(chromosomes[0]);
           #print([i.evaluation_criteria for i in chromosomes])

    else:
       best = copy.deepcopy(chromosomes[0]);       
      
    return chromosomes,best
              
def selection(chromosomes,selection_rate):

    # how many chromosomes do we keep
    n_keep = int(np.round(selection_rate*len(chromosomes)))
   
    # select
    chromosomes = chromosomes[0:n_keep]
    
    return chromosomes,n_keep
    
def pairing(chromosomes,n_chromosomes,n_keep):

    parents = []
    
    for i in range(int(np.ceil((n_chromosomes-n_keep)/2.))):
        
        potential = chromosomes[:]
        first,potential = tournament_selection(potential)       
        second,_ = tournament_selection(potential)
        parents.append([first,second])

    return parents
    
def tournament_selection(potential):

    # first parent candidate
    first =  random.choice(potential)
    potential.remove(first)
    
    # second parent candidate    
    second =  random.choice(potential)
    potential.remove(second)  
    # select the chromosoom with the largest peformance
    if (first.evaluation_criteria < second.evaluation_criteria):
        potential.append(first)
        return second,potential
    else:
        potential.append(second)
        return first,potential

def mating(chromosomes,parents,crossover_rate,mode,n_chromosomes,n):
        
    for i in range(len(parents)):
       
        # get parameters parents
        parameters_parent1 = parents[i][0].parameters
        parameters_parent2 = parents[i][1].parameters
        
        if np.random.rand()<crossover_rate:
            
            # mating binary
            if mode == "binary":
                chromosomes,n = mating_binary(chromosomes,parameters_parent1,parameters_parent2,n)
    
            # matng continues
            if mode == "continuous":
                chromosomes,n = mating_continues(chromosomes,parameters_parent1,parameters_parent2,n)
        
        else:
            
            # copy parents
            new = Chromosoom(parameters_parent1);new.setID(n);n+=1
            chromosomes.append(new)
     
            new = Chromosoom(parameters_parent2);new.setID(n);n+=1
            chromosomes.append(new)
    
    # keep size population constant
    chromosomes[0:n_chromosomes]
    
    return chromosomes,n
    
def mating_binary(chromosomes,parameters_parent1,parameters_parent2,n):
    
    # choose beta (go outside range)
    l = len(parameters_parent1)
    beta = int(np.round(np.random.rand(1)*l))
    
    # copy parameters for offspring 1
    parameters_offspring = parameters_parent1.copy()
    parameters_offspring["sample"].iloc[beta:] =parameters_parent2["sample"].iloc[beta:]
    new = Chromosoom(parameters_offspring);new.setID(n);n+=1
    chromosomes.append(new);
    
    # copy parameters for offspring 2
    parameters_offspring = parameters_parent1.copy()
    parameters_offspring["sample"].iloc[:-(l-beta)] =parameters_parent2["sample"].iloc[:-(l-beta)]
    new = Chromosoom(parameters_offspring);new.setID(n);n+=1
    chromosomes.append(new);
  
    return chromosomes,n

def mating_continues(chromosomes,parameter_parent1,parameter_parent2,n):
    
    cond=True
    beta = np.random.rand(len(parameters_parent1))*1.2-0.1
        
    for i in range(2):
        
        while cond==True:
            
            # simulate parameters offspring
            parameters_offspring = parameters_parent1.copy()
            parameters_offspring["sample"] = beta*(parameters_parent1["sample"]-parameters_parent2["sample"])+parameters_parent1["sample"]

            # apply boundary conditions
            cond = check_boundary_conditions(parameters_offspring)

        new = Chromosoom(parameters_offspring);new.setID(n);n+=1
        chromosomes.append(new);
  
    return chromosomes,n
    
def mutation(chromosomes,mutation_rate,mode,n):
    
    n_mutation = int(np.round(mutation_rate*len(chromosomes[0].parameters)*len(chromosomes)))
    
    for i in range(n_mutation):

        # choose a chromosoom
        chromosoom = random.choice(chromosomes)
        # get parameters and randoimport shutillmize one element
        init_parameters = chromosoom.parameters.copy()
        # remove from the pool
        chromosomes.remove(chromosoom)            
                 
        if mode == "continuous":
            
            mutated,n =  mutation_continuous(init_parameters,n)
        
        if mode == "binary":
            
            mutated,n = mutation_binary(init_parameters,n)
            
        chromosomes.append(mutated)

    return chromosomes,n

def mutation_continuous(init_parameters,n):

    # Iterate untill condition is satisfied
    cond=True
    
    while cond==True:

        # mutate
        parameters =init_parameters.copy()
        parameters["sample"][random.randrange(0,len(parameters))] = np.random.rand()
        # check boundary conditions
        cond = check_boundary_conditions(parameters)

    # make new chromosoomm and add to chromosomes
    mutated = Chromosoom(parameters)
    mutated.setID(n);n+=1
    
    return mutated,n
    
def mutation_binary(init_parameters,n):

    # Iterate untill condition is satisfied

    # mutate
    parameters =init_parameters.copy()
    r = random.randrange(0,len(parameters))
    parameters["sample"][r] = 1. if  parameters["sample"][r] == 0. else 0.

    # make new chromosoomm and add to chromosomes
    mutated = Chromosoom(parameters)
    mutated.setID(n);n+=1
    
    return mutated,n
       
def print_chromosooms(chromosomes,res,run):

    # performance
    perf_keys = list(chromosomes[0].performance.keys());
    # columns    
    columns = list(chromosomes[0].parameters["parameter"]) + ["ID"] + perf_keys
    
    # data
    data = np.zeros([len(chromosomes),len(columns)])
    
    for i in range(len(chromosomes)):
        
        data[i,:-(1+len(perf_keys))] = np.array(chromosomes[i].parameters["sample"]).ravel()
        data[i,-(1+len(perf_keys))]  = chromosomes[i].ID
        data[i,-len(perf_keys):]  = np.array([chromosomes[i].performance[perf_keys[j]] for j in range(len(perf_keys))])
    
    dataframe = pd.DataFrame(data=data,columns=columns)
    #pd.DataFrame(data=data,columns=columns).to_csv(os.path.join(res,str(run)+".csv"))
    dataframe.to_csv(os.path.join(res,str(run)+".csv"))
    
    #f = open(os.path.join(res,"std.txt"),"a")
    #f.write(str(np.std(dataframe["mae"]))+"\n")
    #f.close()

    #f = open(os.path.join(res,"mean.txt"),"a")
    #f.write(str(np.mean(dataframe["mae"]))+"\n")
    #f.close()
    
def multithread(chromosomes,model,model_input,boundaries,objective_function,nan_value,res,ncpu=-1,full_output=False):
    
    # import package
    import multiprocessing
    
    # start timer
    nos = np.sum([1. if np.isnan(chromosomes[i].evaluation_criteria) else 0. for i in range(len(chromosomes))])
    t = Runtime(nos,10)
    # get number of processors
    if ncpu==-1:
        ncpu=multiprocessing.cpu_count()



    
    f = open(os.path.join(res,"cpu.txt"),"a") 
    f.write("Have to calculate "+str(nos)+" instances \n")
    f.write("Starting pool with %s processes" %ncpu+"\n")
    
    if ncpu!=1:
    
        # pool processors
        pool=multiprocessing.Pool(ncpu)
        # make list of jobs
        jobs=[0.]*len(chromosomes)
        for i in range(len(chromosomes)):
            if np.isnan(chromosomes[i].evaluation_criteria):
                jobs[i] = pool.apply_async(eval(model),(model_input,boundaries,chromosomes[i],nan_value,res,full_output))  
        pool.close()
        
        for i in range(len(chromosomes)):
            if np.isnan(chromosomes[i].evaluation_criteria):  
                performance = jobs[i].get()
                chromosomes[i].setPerformance(performance)
                chromosomes[i].setEvaluationcriteria(performance[objective_function])
        pool.join()
  
    else:
        
        for i in range(len(chromosomes)):
             if np.isnan(chromosomes[i].evaluation_criteria):
                 performance= eval(model+"(model_input,boundaries,chromosomes[i],nan_value,res,full_output)")
                 chromosomes[i].setPerformance(performance)
                 chromosomes[i].setEvaluationcriteria(performance[objective_function])
                         
    dt = t.close()
    f.write("End pool, total runtime is of iteration is: "+str(dt)+"\n")
    f.write("------\n")
        
    return chromosomes        

def create_dir(res,L):
    
    for i in range(len(L)):
        if not os.path.exists(os.path.join(res,L[i])):
            os.makedirs(os.path.join(res,L[i]))      
    
    
