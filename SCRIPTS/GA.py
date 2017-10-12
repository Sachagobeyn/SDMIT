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
from copy import deepcopy
import random
from script import *
import shutil
import os
import pandas as pd
#pd.set_option('chained_assignment',None)
from runtime import Runtime
import multiprocessing
from collections import OrderedDict

minimize = ["AIC","SSE","BIC"]

def GA(model,model_inputs,boundaries,options,res,full_output=False):
    
    print("="*21)
    print("Start GA optimisation")
    print("="*21)
    
    " Load hyper parameters of evolutionary algorithm"
    n_chromosomes,selection_rate,mutation_rate,crossover_rate,criteria,duplicates,adaptive,k1,k2,k3,k4 = load_hyper_parameters(options)

    " Set environment EA, type of optimisation problem and objective function"
    mode,multi_objective,objective_function,nan_value,ncpu = set_environment_EA(options,res)
    
    " Write initial logs"
    write_logs(n_chromosomes,selection_rate,mutation_rate,crossover_rate,mode,nan_value,adaptive,k1,k2,k3,k4,multi_objective,duplicates,res)
    
    " Initialize chromosomes and population"
    population,run,ID = initialize_population(model,model_inputs,boundaries,objective_function,n_chromosomes,mode,multi_objective,res,ncpu,nan_value,full_output=full_output)

    " Start timer"
    t = Runtime(n_chromosomes,100)
       
    " Initiate vector containing best solutions for every generation"
    best_criteria = []
    
    "Start"
    cond  = True      
    while cond:
        
        # If preferable (not tested with latest version) remove duplicates + introduce random chromosomes
        #if duplicates==False:
        #    population.chromosomes,ID = proces_duplicates(population.chromosomes,'fitness',boundaries,n_chromosomes,ID) 
        
        # Check the conditions for running GA (i.e. while cond)
        run+=1
        cond = False if run >= options["maximum_runs"] else True    
        
        if (run>criteria):

            cond = False if np.sum(np.array(best_criteria[-(int(criteria)):])==best_criteria[-1])==criteria else True
   
        # If stopping criteria are not met, select, pair, mate and mutate population
        if cond==True:
            
            # Open runtime
            t = Runtime(n_chromosomes,100)

            # If preferable (not tested with latest version) remove duplicates + introduce random chromosomes            
            #if duplicates==False:
            #    population.chromosomes,ID = proces_duplicates(population.chromosomes,'fitness',boundaries,n_chromosomes,ID) 
            
            # Step 1: Selection 
            parents = selection(population.chromosomes,selection_rate,multi_objective)
                         
            # Step 2: Cross-over
            parents,offspring,ID = crossover(parents,crossover_rate,mode,n_chromosomes,ID,multi_objective,adaptive=adaptive,fmean=population.meanPrintCrit(nan_value),fmax=population.best.fitness,k1=k1)
            
            # Step INTER: did testing with adaptive parameters, however, removed from final code
#            if adaptive==True:
#                population.setPopulation(chromosomes)
#                population.calculateChromosomes(ncpu,os.path.join(res,"model_runs"),full_output,mode,nan_value,multi_objective)
            
            # Step 3: Mutation
            offspring,ID = mutation(offspring,mutation_rate,mode,ID,adaptive=adaptive,fmean=population.meanPrintCrit(nan_value),fmax=population.best.fitness,k2=k2)
            
            # Step 4: Calculate new chromosomes
            population.setPopulation(parents + offspring)
            population.calculateChromosomes(ncpu,os.path.join(res,"model_runs"),full_output,mode,nan_value,multi_objective)

            # Step 5: Evaluate and save
            chromosomes,best =  evaluate_chromosomes(population.chromosomes,population.best)
            population.setPopulation(chromosomes)
            population.setBest(best)
            
            # Step 6: Determine fronts + density + resize population
            if multi_objective==True:
                population.FNDS()
                population.CDA()  
                population.resize()
            
            # Step 7: Dynamically adapt boundaries for optimisation (losen boundaries)
            "Only for parameter estimation"
            if mode!="binary":
                population.adjustBoundaries(criteria,mode,population.best,run,multi_objective,os.path.join(res,"boundaries"),nan_value)
            
            # Step 8: Print population and print message for generation
            if cond==True:
                population.printPopulation(os.path.join(res,"optimisation_summary"),run,mode,multi_objective)
#            population.printPopulation(os.path.join(res,"optimisation_summary"),run)
            best_criteria.append(population.best.fitness)
            dt = t.close(print_value=False)    
            print_message(population,run,dt,multi_objective,res)
    
        
    print("="*19)
    print("End GA optimisation")
    print("="*19)
       
    return population.best.performance,population.chromosomes

def load_hyper_parameters(options):
    
    n_chromosomes  = int(options["number_of_chromosomes"]);
    selection_rate = float(options["selection_rate"])
    mutation_rate  = float(options["mutation_rate"])
    crossover_rate = float(options["crossover_rate"])
    criteria = options["stop_criteria"]
    duplicates = options["duplicates"]

    " TO DO: implement conditions"
    if n_chromosomes<8:
        print(" Population too small (<8), changing to standard (32)")
        n_chromosomes = 32
        
    adaptive = options["adaptive"]
    k1 = float(options["k1"])
    k2 = float(options["k2"])
    k3 = float(options["k3"])
    k4 = float(options["k4"]) 

    adaptive = False
    k1 = np.nan
    k2 = np.nan
    k3 = np.nan 
    k4 = np.nan
    
    return n_chromosomes,selection_rate,mutation_rate,crossover_rate,criteria,duplicates,adaptive,k1,k2,k3,k4

def set_environment_EA(options,res):
 
    " create directory "
    l = ["optimisation_summary","model_runs","boundaries"]
    for i in l:
        if os.path.exists(os.path.join(res,i)):
            shutil.rmtree(os.path.join(res,i))    
    create_dir(res,l)
    
    " Type of problem"
    mode = options["mode"]
     
    " Objective function"
    multi_objective = options["multi_objective"]
    objective_function = options["objective_function"]
    
    " Nan values un print"
    nan_value = options["nan_value"]

    " Multiprocessing"

    if "ncpu" in options:
        ncpu = int(options["ncpu"])
    else:
        ncpu = -1

    return mode,multi_objective,objective_function,nan_value,ncpu
    
        
class Population():

    def __init__(self,chromosomes):
        
        self.chromosomes = chromosomes
        self.n = len(self.chromosomes)
        self.best = np.nan
        
    def addchromosome(self,chromosome):
        
        self.chromosomes.append(chromosome)
        self.n = len(self.chromosomes)

    def deletechromosome(self,chomosoom):

        self.chromosomes.remove(chromosome)      
    
    def setPopulation(self,chromosomes):
        
        self.chromosomes = chromosomes
    
    def setBest(self,best):
        
        self.best = best
        
    def setEvaluation(self,model,model_input,boundaries,objective_function):
        
        self.model = model
        self.model_input = model_input
        self.boundaries = boundaries
        self.objective_function = objective_function
       
    def adjustBoundaries(self,criteria,mode,best,run,multi_objective,resmap,nan_value):

        for i in self.chromosomes:
        
            if multi_objective==False:
            
                if i.fitness==best.fitness:
                    
                    i.mapParameters(np.nan)
                    # only for boundaries of continuous variables
                    cond = self.boundaries["type"]=="continuous"
                    if np.sum(cond)>0:
                        parameters_i = i.parameters[cond]
                        boundaries = self.boundaries[cond]
                        # extent b4   
                        cond_b4 = (parameters_i["a4"]>boundaries["b4"].astype(float)) & (~np.isnan(parameters_i["a4"]))
                        boundaries.loc[cond_b4,"b4"] = parameters_i["a4"][cond_b4].values
                        # extent b3
                        cond_b3 = (parameters_i["a3"]<boundaries["b3"].astype(float)) & (~np.isnan(parameters_i["a3"]))
                        boundaries.loc[cond_b3,"b3"] = parameters_i["a3"][cond_b3].values
                        # extent b3
                        cond_b2 = (parameters_i["a2"]>boundaries["b2"].astype(float)) & (~np.isnan(parameters_i["a2"]))
                        boundaries.loc[cond_b2,"b2"] = parameters_i["a2"][cond_b2].values                        
                        # extent b1
                        cond_b1 = (parameters_i["a1"]<boundaries["b1"].astype(float)) & (~np.isnan(parameters_i["a1"]))
                        boundaries.loc[cond_b1,"b1"] = parameters_i["a1"][cond_b1].values
                        
                        self.boundaries[cond] = boundaries                  
            
            else:
            
                if i.rank==1:
                    
                    i.mapParameters(np.nan)
                    # only for boundaries of continuous variables
                    cond = self.boundaries["type"]=="continuous"
                    if np.sum(cond)>0:
                        parameters_i = i.parameters[cond]
                        boundaries = self.boundaries[cond]
                        # extent b3                
                        cond_b4 = (parameters_i["a4"]>boundaries["b4"].astype(float)) & (~np.isnan(parameters_i["a4"]))
                        boundaries.loc[cond_b4,"b4"] = parameters_i["a4"][cond_b4].values
                        # extent b0
                        cond_b1 = (parameters_i["a1"]<boundaries["b1"].astype(float)) & (~np.isnan(parameters_i["a1"]))
                        boundaries.loc[cond_b1,"b1"] = parameters_i["a1"][cond_b1].values
                        
                        self.boundaries[cond] = boundaries    
        
        self.boundaries.to_csv(os.path.join(resmap,"boundaries_"+str(run)+".csv"))
        
    def remove_duplicates(self,attribute):
        
        seen = list()
        unique = []
        for obj in self.chromosomes:
            if eval('obj.'+str(attribute)) not in seen:
                unique.append(obj)
                seen.append(eval('obj.'+str(attribute)))
        return unique
        
    def proces_duplicates(self,attributes,boundaries,n_chromosomes,ID):
        
        self.remove_duplicates(attributes)
        
        for i in range(n_chromosomes-len(self.chromosomes)):
            
            new,ID = initiate_chromosome(boundaries,ID)
            self.chromosomes.append(new)
                
        return ID
        
    def meanPrintCrit(self,nan_value):
        
        return np.nanmean([i.printCrit for i in self.chromosomes if i.printCrit !=-nan_value])
        
    def stdPrintCrit(self,nan_value):
        
        return np.nanstd([i.printCrit  for i in self.chromosomes if i.printCrit !=-nan_value])

    #fast-non-dominated-sort
    def FNDS(self):
    
        #Deb, K., Pratap, A., Agarwal, S., Meyarivan, T., 2002. A fast and elitist 
        #multiobjective genetic algorithm: NSGA-II. IEEE Trans. Evol. Comput. 6, 182–197.
        self.F = [[]]
        S = {p.ID:[] for p in self.chromosomes}
        n = {p.ID:0 for p in self.chromosomes}
        
        # first front    
        for p in range(len(self.chromosomes)):
                    
            for q in range(len(self.chromosomes)):
                
                if dominates(self.chromosomes[p].ofs,self.chromosomes[q].ofs)==True:
                    
                    S[self.chromosomes[p].ID]  += [self.chromosomes[q]]
                    
                elif dominates(self.chromosomes[q].ofs,self.chromosomes[p].ofs)==True:
                    
                    n[self.chromosomes[p].ID] += 1
                    
            if n[self.chromosomes[p].ID]==0:
                self.chromosomes[p].rank = 1
                self.F[0] += [self.chromosomes[p]]
               
                
        i = 1
        # other fronts    
        
        while len(self.F[i-1])>0:
            
            Q = []
           
            for p in self.F[i-1]:
                
                for q in S[p.ID]:
                    
                    n[q.ID] = n[q.ID]-1  
                    
                    if n[q.ID]==0:
                        
                        self.chromosomes[[chromosome.ID==q.ID for chromosome in self.chromosomes].index(True)].rank = i+1
                        q.rank = i+1
                        Q += [q]
            i = i+1
            self.F += [Q]
        
        self.F = self.F[:i-1]

    # crowded distance measure
    def CDA(self,error=10**-10):
        #Deb, K., Pratap, A., Agarwal, S., Meyarivan, T., 2002. A fast and elitist 
        #multiobjective genetic algorithm: NSGA-II. IEEE Trans. Evol. Comput. 6, 182–197.
        
        for I in self.F:
            
            D = {i.ID:0 for i in I}
            
            # for number of objectives
            for m in range(len(I[0].ofs)):
                
                of_m = {i.ID:i.ofs[m] for i in I}
                
                # sort 
                of_m = OrderedDict(sorted(of_m.items(), key=lambda x: x[1]))
                keys = list(of_m.keys())
                
                D[keys[0]] = 10**10
                D[keys[-1]] = 10**10
                
                for i in range(0,len(keys)-1,1):
            
                    # "+ error": compensate if np.min and np.max in fraction are 0.
                    D[keys[i]] = D[keys[i]] + (of_m[keys[i+1]]-of_m[keys[i-1]])/(np.max(list(of_m.values()))-np.min(list(of_m.values()))+error)
            
            for i in I:
                
                self.chromosomes[[chromosome.ID==i.ID for chromosome in self.chromosomes].index(True)].distance = D[i.ID]
      
    def resize(self):

        # step 1: sort population according to rank
        self.chromosomes.sort(key=operator.attrgetter('distance')) 
        self.chromosomes.sort(key=operator.attrgetter('rank'))
                  
        # step 2: select first n chromosomes
        self.chromosomes = self.chromosomes[0:self.n]
     
    def calculateChromosomes(self,ncpu,res,full_output,mode,nan_value,multi_objective,final_run=False):
        
        "Multi-/singlethread"
        # make list of jobs
        if ncpu!=1:
 
            if ncpu==-1:
                
                ncpu=multiprocessing.cpu_count()
            
            self.multiProcess(ncpu,res,full_output,mode,nan_value,multi_objective,final_run=final_run)
    
        else:       
            
            self.singleProcess(res,full_output,mode,nan_value,multi_objective,final_run=final_run)
                        
    def multiProcess(self,ncpu,res,full_output,mode,nan_value,multi_objective=False,final_run=False):
        
        # Make pool
        pool=multiprocessing.Pool(ncpu)
        # Make jobs
        jobs=[0.]*len(self.chromosomes)
        
        for i in range(len(self.chromosomes)):
            if np.isnan(self.chromosomes[i].fitness):
                jobs[i] = pool.apply_async(eval(self.model),(self.model_input,self.boundaries,self.chromosomes[i],mode,nan_value,res,full_output,final_run))
        pool.close()
        
        for i in range(len(self.chromosomes)):
            if np.isnan(self.chromosomes[i].fitness):  
                performance,solution = jobs[i].get()
                self.chromosomes[i].setPerformance(performance)
                self.chromosomes[i].setSolution(solution)
                self.chromosomes[i].setPrintCrit(performance["TSS"])
                if type(self.objective_function)==list:
                    self.chromosomes[i].setOFs([performance[j] for j in self.objective_function])
                else:
                    self.chromosomes[i].setFitness(-performance[self.objective_function] if self.objective_function in minimize else performance[self.objective_function])
        pool.join()
        
    def singleProcess(self,res,full_output,mode,nan_value,multi_objective=False,final_run=False):
 
         for i in range(len(self.chromosomes)):
             if np.isnan(self.chromosomes[i].fitness):
                 performance,solution = eval(self.model+"(self.model_input,self.boundaries,self.chromosomes[i],mode,nan_value,res,full_output,final_run)")
                 self.chromosomes[i].setPerformance(performance)
                 self.chromosomes[i].setSolution(solution)
                 self.chromosomes[i].setPrintCrit(performance["TSS"])
                 if type(self.objective_function)==list:
                     self.chromosomes[i].setOFs([performance[j] for j in self.objective_function])
                 else:
                    self.chromosomes[i].setFitness(-performance[self.objective_function] if self.objective_function in minimize else performance[self.objective_function])

         return performance,solution
        
    def printPopulation(self,res,run,mode,multi_objective,initiate=False):

        "print to seperate file"
        dataframe,parameters = self.transformToOutput(mode,multi_objective)
        #pd.DataFrame(data=data,columns=columns).to_csv(os.path.join(res,str(run)+".csv"))
        dataframe.to_csv(os.path.join(res,str(run)+".csv"))
        
        "initiate/append history of algorithm"
        dataframe.loc[:,"generation"] = run
        parameters.loc[:,"generation"] = run
        
        if initiate==True:
            
            dataframe.to_csv(os.path.join(res,"results_iterations.csv"))
            parameters.to_csv(os.path.join(res,"model_parameters.csv"))
            
        else:
        
            dataframe.to_csv(os.path.join(res,"results_iterations.csv"),mode='a',header=False)
            parameters.to_csv(os.path.join(res,"model_parameters.csv"),mode='a',header=False)
           
    def transformToOutput(self,mode,multi_objective):
         
        
        un_var = self.chromosomes[0].parameters["variable"].unique()
        un_var.sort()
        
        # performance
        perf_keys = list(self.chromosomes[0].performance.keys());perf_keys.sort()
        columns = list(un_var)+["ID"]+["mutated"]+perf_keys+["rank"]
        # data
        data = np.full([len(self.chromosomes),len(columns)],None, dtype=np.object)
        
        for i in range(len(self.chromosomes)):
            
            self.chromosomes[i].parameters = self.chromosomes[i].parameters.sort_values("variable")
            if (mode=="variable") or (mode=="continuous"):
                data[i,:-(3+len(perf_keys))] = [0 if (type(j)==int) | (type(j)==np.int64) else " ".join(["%.3f"%k for k in j.returnString()]) for j in self.chromosomes[i].parameters["sample"]]
            elif mode=="binary":
                data[i,:-(3+len(perf_keys))] = np.array(self.chromosomes[i].parameters["sample"]).ravel()  
            data[i,-(3+len(perf_keys))]  = self.chromosomes[i].ID
            data[i,-(2+len(perf_keys))]  = self.chromosomes[i].mutated
            data[i,-(1+len(perf_keys)):-1]  = np.array([self.chromosomes[i].performance[perf_keys[j]] for j in range(len(perf_keys))])
            if multi_objective==True:
                data[i,-1] = self.chromosomes[i].rank

        # get best solution
        self.best.mapParameters(np.nan)
        self.best.parameters.loc[:,"ID"] = deepcopy(self.best.ID)
       
        return pd.DataFrame(data=data,columns=columns)[columns],self.best.parameters
    
class chromosome():
    
    def __init__(self,parameters):
        
        self.parameters = parameters
        self.performance = np.nan
        self.printCrit = np.nan
        self.ID = np.nan
        self.fitness = np.nan
        self.solution = np.nan
        self.mutated= False 
        
        # Information for non-dominated sorting
        self.f = np.nan
        self.ofs = np.nan     
        self.rank = np.nan
        self.distance = np.nan
        self.protect = False
    
    def setProtect(self,protect):
        
        self.protect = protect
        
    def setPerformance(self,performance):
        
        self.performance = performance
    
    def setFitness(self,fitness):
        
        self.fitness = fitness
    
    def setPrintCrit(self,printCrit):
        
        self.printCrit = printCrit
        
    def setOFs(self,ofs):
        
        self.ofs = ofs

    def setID(self,n):
        
        self.ID = n
        
    def setSolution(self,solution):
        
        self.solution = solution
        
    def setParameters(self,parameters):
        
        self.parameters = parameters
 
    def mapParameters(self,nan_value):
        
        index = self.parameters.index

        for i in index:
            
            if type(self.parameters.loc[i,"sample"])!=int:
             
                 a = ["a"+str(j) for j in range(1,len(self.parameters.loc[i,"sample"].returnString())+1,1)]
                 
                 for j in range(len(a)):

                    self.parameters.loc[i,a[j]] = self.parameters.loc[i,"sample"].parameters[j]
                
            else:
                
                a = ["a1","a2","a3","a4"]
                
                for j in a:
                    
                    self.parameters.loc[i,j] = nan_value
                            
    def mutation_operator(self,mutation_rate,mode,ID,adaptive=False,**kwargs):
        
        self.mutated = False
        
        if adaptive==True:
            
            fp = self.fitness
            mutation_rate = kwargs["k2"]*(kwargs["fmax"]-fp)/(kwargs["fmax"]-kwargs["fmean"]) if fp>=kwargs["fmean"] else kwargs["k2"]        
    
        mutation_flag = False
        
        if 0<mutation_rate:
         
            mutation_flag = True            
            self.setID(deepcopy(ID))
            
            if mode == "variable":
                
                self.mutation_variable(mutation_rate)
            
            if mode == "binary":
                
                self.mutation_binary(mutation_rate)
        
            if mode == "continuous":
            
                self.mutation_continuous(mutation_rate)
                
            self.mutated = True
            
        return mutation_flag
        
#    def mutation_continuous(self,ID):
#    
#        print("To be coded")
#        # Iterate untill condition is satisfied
#        cond=True
#        
#        while cond==True:
#    
#            # mutate
#            parameters =init_parameters.copy()
#            parameters["sample"][random.randrange(0,len(parameters))] = np.random.rand()
#            # check boundary conditions
#            cond = check_boundary_conditions(parameters)
#    
#        # make new chromosomem and add to chromosomes
#        ID+=1
#        
#        return parameters,ID
#        return chromosome
        
    def mutation_binary(self,mutation_rate):
        
        # mutate
        index = self.parameters.index
        r = random.randrange(0,len(index))        
        self.parameters.loc[index[r],"sample"] = 1. if self.parameters.loc[index[r],"sample"] == 0. else 0.
        
        # set new properties chromosome
        self.fitness = np.nan
        self.rank  = np.nan
        self.distance = np.nan

    def mutation_variable(self,mutation_rate):
                    
        index = self.parameters.index

        for i in index:

            "mutate the parameters OR mutate the presence/absence of a parameter"
        
            if type(self.parameters.loc[i,"sample"])!=int:
                
                "Switch parameter off with probability equal to the mutation rate"
                if (np.random.uniform()>mutation_rate):
                                            
                    self.parameters.loc[i,"sample"].manipulate(mutation_rate)
                    
                else:
                    
                    self.parameters.loc[i,"sample"] = 0
            else:
                
                if (np.random.uniform()<mutation_rate):
                                    
                    self.parameters.loc[i,"sample"] = substring(self.parameters.loc[i,["b1","b2","b3","b4"]],self.parameters.loc[i,"b1"],self.parameters.loc[i,"b4"],4,self.parameters.loc[i,"type"],initiate=True)


        self.fitness = np.nan
        self.rank = np.nan
        self.distance = np.nan
        
    def mutation_continuous(self,mutation_rate):
    
        "mutate every substring bit with equal probability 'mutation_rate'"
        index = self.parameters.index	

        for i in index:
            
            self.parameters.loc[i,"sample"].manipulate(mutation_rate)

        self.fitness = np.nan
        self.rank = np.nan
        self.distance = np.nan
           
    def deepcopy_parameters(self):

        parameters = deepcopy(self.parameters)
        index = self.parameters.index	
        
        for i in index:
        
            if type(self.parameters["sample"].loc[i])!=int:
                
                parameters.loc[i,"sample"] = deepcopy(self.parameters.loc[i,"sample"])

        return parameters             
                 
    def calculatechromosome(self,model,model_inputs,boundaries,objective_function,nan_value,res,full_output=False,final_run=True):
 
         self.performance,self.solution = eval(model+"(model_inputs,boundaries,self,nan_value,res,full_output,final_run)")
         self.setFitness(-self.performance[objective_function])  

    def writeparameterschromosome(self,res):
        
        self.parameters["sample"] = [self.parameters["sample"].iloc[j] if type(self.parameters["sample"].iloc[j])==int else np.array2string(self.parameters["sample"].iloc[j].returnString()) for j in range(len(self.parameters))]
        self.parameters.to_csv(res)
        
#class bitstring():
#
#    # note: Grey encoding, possible to encode to grey and int    
#    def __init__(self,precision=4):
#        
##        self.precision = int(precision)        
##        self.gstring = np.array([str(random.randint(0,1)) for i in range(self.precision)],dtype = int)
##        self.values = np.ones([self.precision])
##        for i in range(self.precision):
##            self.values[i] = self.values[i-1]*2 if i!=0 else 1
##        self.computeInt()
#        
#        # define precision
#        self.precision =  int(precision)
#        # define random int
#        self.values = np.ones([self.precision])
#        for i in range(self.precision):
#            self.values[i] = self.values[i-1]*2 if i!=0 else 1
#        self.int = int(np.sum([random.randint(0,1)*self.values[i] for i in range(self.precision)]))
#        self.computeString()
#         
#    def int2bin(self):
#
#        n = deepcopy(self.int)
#        ind = 1
#        
#        if n:
#            self.bstring = [0]*self.precision
#            while n:
#                n,remainder = divmod(n, 2)
#                self.bstring[-ind] = int(remainder)
#                ind +=1
#            
#        else: 
#            self.bstring = [0]*self.precision
#             
#    def bin2int(self):
#
#        i = 0
#        for bit in self.bstring:
#            i = i * 2 + bit
#        self.int = i
#
#    def bin2grey(self):
#        
#        self.gstring = self.bstring[:1] + [1 if self.bstring[i-1]!=self.bstring[i] else 0 for i in range(1,len(self.bstring))]
#        
#    def grey2bin(self):
#        
#        bstring = [0]*self.precision
#        bstring[0] = self.gstring[0]
#        
#        for i in range(1,len(self.bstring)):
#            if self.gstring[i]==0:
#                bstring[i] = bstring[i-1]
#            else:
#                bstring[i] = abs(bstring[i-1]-1)
#        
#        self.bstring = bstring
#        
#    def computeInt(self):
#        
#        self.grey2bin()
#        self.bin2int()
#
#    def computeString(self):
#
#        #convert to bstring
#        self.int2bin()
#        #convert to gstring
#        self.bin2grey()
#        
#    def manipulate(self,chance):
##        
#        cond = np.array([np.random.uniform()<chance for i in range(len(self.gstring))])
#        newstring = np.array(self.gstring)
#        newstring[cond] = np.abs(newstring[cond]-1)
#        self.gstring = list(newstring)
#        self.computeInt()
###     
#    def setInt(self,int_):
#
#        self.int = int(int_)
#        self.computeString()
#        
##    def setString(self,bstring):
##
###        self.bstring =  bstring
##
##        
#    def returnInt(self):
#        
#        self.computeInt()
#        
#        return self.int 
#        
##    def returnString(self):
###        
###        self.computegString()
###        
###        return self.string 
#    

class substring():
    
    def __init__(self,boundaries,low,high,n,ctype,initiate=False):
        
        self.boundaries = boundaries
        self.type = ctype
        if (initiate==True) & (ctype=="continuous"):
            self.parameters = [np.random.uniform(boundaries[0],boundaries[1]) for i in range(2)]+[np.random.uniform(boundaries[2],boundaries[3]) for i in range(2)]
        if (initiate==True) & (ctype=="categorical"):
            self.parameters = np.array([np.random.uniform(low,high) for i in range(len(boundaries[~boundaries.isnull()]))])
        self.ctype = ctype
        if ctype=="continuous":
            self.parameters = np.sort(self.parameters)#np.array([0.99*self.parameters[i+1] if self.parameters[i+1]<self.parameters[i] else self.parameters[i] for i in range(3)]+[self.parameters[3]])
        self.low = low
        self.high = high
                
    def returnString(self):
        
        return self.parameters
    
    def setString(self,parameters):
        
        self.parameters = parameters
        
    def manipulate(self,chance):
#        
        if self.type=="continuous":
            
            if np.random.uniform()<chance:
                #print("low:%.2f,a2:%.2f,a3:%.2f,high:%.2f"%((self.low,self.parameters[0],self.parameters[1],self.high)))               
                self.parameters[0] = np.random.uniform(self.low*0.5 if self.low>0 else self.low*1.5,self.parameters[1])
                #print("low:%.2f,a2:%.2f,a3:%.2f,high:%.2f"%((self.low,self.parameters[0],self.parameters[1],self.high)))               
            if np.random.uniform()<chance:
                #print("low:%.2f,a2:%.2f,a3:%.2f,high:%.2f"%((self.low,self.parameters[0],self.parameters[1],self.high)))               
                self.parameters[1] = np.random.uniform(self.parameters[0],self.parameters[2])
              #print("low:%.2f,a2:%.2f,a3:%.2f,high:%.2f"%((self.low,self.parameters[0],self.parameters[1],self.high)))               
            if np.random.uniform()<chance:
                #print("low:%.2f,a2:%.2f,a3:%.2f,high:%.2f"%((self.low,self.parameters[0],self.parameters[1],self.high)))               
                self.parameters[2] = np.random.uniform(self.parameters[1],self.parameters[3])
                #print("low:%.2f,a2:%.2f,a3:%.2f,high:%.2f"%((self.low,self.parameters[0],self.parameters[1],self.high)))  
            if np.random.uniform()<chance:
                #print("low:%.2f,a2:%.2f,a3:%.2f,high:%.2f"%((self.low,self.parameters[0],self.parameters[1],self.high)))               
                self.parameters[3] = np.random.uniform(self.parameters[2],self.high*1.5 if self.high>0 else self.high*0.5)
                #print("low:%.2f,a2:%.2f,a3:%.2f,high:%.2f"%((self.low,self.parameters[0],self.parameters[1],self.high)))                  
        else:
            
            for i in range(len(self.parameters)):
                    
                    if np.random.uniform()<chance:
                        
                        self.parameters[i] = np.random.uniform(self.low,self.high)
    
    def fixBoundary(self):
        
        self.parameters = [self.low if self.parameters[i]<self.low else self.parameters[i] for i in range(len(self.parameters))] 
        self.parameters = [self.high if self.parameters[i]>self.high else self.parameters[i] for i in range(len(self.parameters))] 
        if self.type=="continuous":
            self.parameters = np.sort(self.parameters)
        else:
            self.parameters = np.array(self.parameters)
        
def boundaryStringBit(bitString1,bitString2):
    
    "Compare two bitstrings and make sure bitString1 is always <= bitString2"
    bitInt1 = bitString1.returnInt()
    bitInt2 = bitString2.returnInt()
    
    if (bitInt2<bitInt1):
        
        "either we down bitInt1 or rease bitInt2"
        if np.random.uniform()<0.5:
            bitInt1 =  bitInt2 - 1 if bitInt2!=0 else 0
            bitString1.setInt(bitInt1)
        else:
            bitInt2 =  bitInt1 + 1 if bitInt1!=np.sum(bitString1.values) else np.sum(bitString1.values)
            bitString2.setInt(bitInt2)

    return [bitString1,bitString2]
            
            
def write_logs(n_chromosomes,selection_rate,mutation_rate,crossover_rate,mode,nan_value,adaptive,k1,k2,k3,k4,multi_objective,duplicates,res):
    
    print(" Running GA with")
    print("\t "+str(n_chromosomes)+" chromosomes")
    print("\t a selection rate of %.2f"%(selection_rate))
    if adaptive==True:
        print("\t adaptive crossover & mutation operators \n"+"\t"*2 +"k1 = %.2f  \n"%k1+"\t"*2 +"k2 = %.2f \n"%k2 +"\t"*2 +"k3 = %.2f \n"%k3+"\t"*2 +"k4 = %.2f)"%k4) 
    else:
        print("\t crossover & mutation operators \n"+"\t"*2 +"p_c = %.2f \n"%crossover_rate+"\t"*2 + "p_m = %.2f"%mutation_rate)
    
    #print(" Local optimisation: %s"%str(local_optimisation))
    print(" Duplicates allowed: %s"%str(duplicates))
    print(" Multi-objective optimisation %s"%str(multi_objective))
    print(" Writing to .\%s"%res)
    print("="*19)
                    
def initialize_chromosomes(boundaries,n_chromosomes,mode,ID,nan_value):
    
    chromosomes = [0.]*n_chromosomes
    
    for i in range(n_chromosomes):
    
        # list of lists encoding for MI
        if mode == "variable":
            parameters = sample_parameters_variable(boundaries,nan_value)
        # binary encoding for IVS
        if mode == "binary":
            parameters = sample_parameters_binary(boundaries)
        # continuous encoding for PE
        if mode == "continuous":
            parameters = sample_parameters_continuous(boundaries,nan_value)

        chromosomes[i] = chromosome(parameters)
        chromosomes[i].setID(deepcopy(ID));ID=ID+1

    return chromosomes,ID

def initialize_population(model,model_inputs,boundaries,objective_function,n_chromosomes,mode,multi_objective,res,ncpu,nan_value,full_output=True):
    
    # Set first ID
    ID = 0
    t = Runtime(n_chromosomes,100.)
    # Initialize chromosomes
    chromosomes,ID = initialize_chromosomes(boundaries[["variable","type","low","b1","b2","b3","b4","high","a1","a2","a3","a4"]],n_chromosomes,mode,ID,nan_value)
    # Initialize population
    population = Population(chromosomes)
    # Initialize details for calculating performance
    population.setEvaluation(model,model_inputs,boundaries,objective_function)
    # Evaluate chromosomes
    population.calculateChromosomes(ncpu,os.path.join(res,"model_runs"),full_output,mode,nan_value,multi_objective)
    population.setPopulation(chromosomes)
    chromosomes,best = evaluate_chromosomes(population.chromosomes,population.best)
    population.setBest(best)
    # Non dominated sorting
    if multi_objective==True:
        population.FNDS()
        population.CDA()
    # Cloe initialisation and print files
    dt = t.close(print_value=False)    
    # print first batch
    run = 0 
    print_message(population,run,dt,multi_objective,res)
    population.printPopulation(os.path.join(res,"optimisation_summary"),run,mode,multi_objective,initiate=True)
    
    return population,run,ID


def proces_duplicates(chromosomes,attributes,boundaries,n_chromosomes,n):
    
    chromosomes = remove_duplicates(chromosomes,attributes)
    
    for i in range(n_chromosomes-len(chromosomes)):
        
        new,n = initiate_chromosome(boundaries,n)
        chromosomes.append(new)
        
    return chromosomes,n

def remove_duplicates(objects,attribute):
    
    seen = []
    unique = []
    for obj in objects:
        print(eval('obj.'+str(attribute)))
        if eval('obj.'+str(attribute)) not in seen:
            unique.append(obj)
            seen.append(eval('obj.'+str(attribute)))
    return unique
        
def sample_parameters_variable(boundaries,nan_value):
    
    parameters = boundaries.copy(deep=True)
    parameters["sample"] = np.random.randint(2.,size=len(parameters))
    
    for i in ["a1","a2","a3","a4"]:
        
        parameters.loc[:,i] = nan_value
        
    "add bits for accute"
    index = parameters.index
    
    for i in index:
    
        if parameters.loc[i,"sample"]==1:
            
            parameters.loc[i,"sample"] = substring(parameters.loc[i,["b1","b2","b3","b4"]],parameters.loc[i,"b1"],parameters.loc[i,"b4"],4,parameters.loc[i,"type"],initiate=True)
                
    return parameters
    
def sample_parameters_binary(boundaries):
    
    parameters = boundaries.copy(deep=True)
    parameters["sample"] = np.random.randint(2.,size=len(parameters))
    
    return parameters

    
def sample_parameters_continuous(boundaries,nan_value):
    
    parameters = boundaries.copy(deep=True)
    parameters["sample"] = 1.
    
    for i in ["a1","a2","a3","a4"]:
        
        parameters[i] = nan_value
        
    "add bits for accute"
    index = parameters.index
    
    for i in index:
    
        parameters.loc[i,"sample"] = substring(parameters.loc[i,["b1","b2","b3","b4"]],parameters.loc[i,"b1"],parameters.loc[i,"b4"],4,parameters.loc[i,"type"],initiate=True)
                
    return parameters


def check_boundary_conditions(parameters):
    
    # fix out of boundar
    parameters.loc[parameters["sample"]<0,"sample"] = np.random.rand()
    parameters.loc[parameters["sample"]>1,"sample"] = np.random.rand()
    
    # check if conditions apply
    index = parameters.index
    for i in index:
        
        if parameters.loc[i,"cond"] !="":
            
            cond = parameters.loc[i,"cond"].split(" ")
            
            if cond[0]=="<":

                if parameters.loc[i,"cond"] > parameters.loc[parameters["parameter"]==cond[1],"sample"]:
                    
                    return True
                
            else:
                
                if parameters.loc[i,"cond"] < parameters.loc[parameters["parameter"]==cond[1],"sample"]:
        
                    return True
    
        else:
            
            return False

def initiate_chromosome(boundaries,ID):
    
    # generate n_chromosomes-n_keep random samples
    parameters = sample_parameters_binary(boundaries)
    new = chromosome(parameters)
    new.setID(ID);ID=ID+1
    
    return new,ID

def evaluate_chromosomes(chromosomes,best):
    
    # sort chromosomes on performance
    chromosomes.sort(key=operator.attrgetter('fitness'),reverse=True)

    # elitism
    if type(best)!=float:
        
        # if there is a chromosome with a better peformanca: replace
        if chromosomes[0].fitness < best.fitness:
            chromosomes.append(best)
            # sort chromosomes on performance
            chromosomes.sort(key=operator.attrgetter('fitness'),reverse=True)
            #print([i.fitness for i in chromosomes])
        else:
           # keep best solution
           best = deepcopy(chromosomes[0]);
           
        chromosomes[0].setProtect(True)
           #print([i.fitness for i in chromosomes])

    else:
       best = deepcopy(chromosomes[0]);       
      
    return chromosomes,best

def dominates(f1,f2):    
# Definition: A solution v(f1) is said to dominate another solution v(f2) if 
#  1. The solution v(f1) is no worse than v(x2) in all objectives
#   AND
#  2. The solution v(f1) is strictly better than v(f2) in at least one objective
#  REF: Deb, K., 2015. Multi-Objective Evolutionary Algorithms. In: Kacprzyk, J., Pedrycz, W. (Eds.), Springer Handbook of Computational Intelligence. Springer Science+BusinessMedia, Berlin, pp. 995–1015.
    # v(f1) is no worse than all objectives
    if np.sum([1 if f1[i] >= f2[i] else 0 for i in range(len(f1))])==len(f1): 
        # v(f1) is strictly better in one 
        if np.sum([1 if f1[i] > f2[i] else 0 for i in range(len(f1))])>0:
            return True
        else:
            return False 
    else:
        return False
            
def selection(chromosomes,selection_rate,multi_objective):

    # how many chromosomes do we keep
    if multi_objective==True:
        selection_rate == 1.
        
    n_keep = int(np.round(selection_rate*len(chromosomes)))
    # select
    #parents = chromosomes[0:n_keep]
    
    # selection based on tournament selection
    parents = []
    
    for i in range(n_keep):
       
        potential = chromosomes[:]
        parent_i,_ = tournament_selection(potential,multi_objective)
        parents.append(deepcopy(parent_i)) #deepcopy because we don't want the parents to have the same adress in the memory since instances now are later seperatly altered (e.g. mutation)   
        
    return parents
    
def pairing(chromosomes,n_chromosomes,n_keep,multi_objective):

    parents = []
    
    for i in range(int(np.ceil((n_chromosomes-n_keep)/2.))):
        
        potential = chromosomes[:]
        first,potential = tournament_selection(potential,multi_objective)       
        second,_ = tournament_selection(potential)
        parents.append([first,second])

    return parents
    
def tournament_selection(potential,multi_objective=False):

    # first parent candidate
    first =  random.choice(potential)
    potential.remove(first)
    
    # second parent candidate    
    second =  random.choice(potential)
    potential.remove(second)  
    
    # select the chromosome based on fitness (single objective)
    if multi_objective!=True:
        if (first.fitness > second.fitness):
            potential.append(second)
            return first,potential
        else:
            potential.append(first)
            return second,potential    
    else:
        #Deb, K., Pratap, A., Agarwal, S., Meyarivan, T., 2002. A fast and elitist 
        #multiobjective genetic algorithm: NSGA-II. IEEE Trans. Evol. Comput. 6, 182–197.
        if (first.rank < second.rank) or ((first.rank==second.rank) and (first.distance > second.distance)):
            potential.append(second)
            return first,potential
        else:
            potential.append(first)
            return second,potential               
       
def crossover(parents,crossover_rate,mode,n_chromosomes,ID,multi_objective,adaptive=False,**kwargs):

    offspring = [];
    
    if multi_objective==True:
        
        cond = len(offspring) < len(parents)
        
    else:
        
        cond = len(offspring)<n_chromosomes-len(parents)
        
    while cond:
        
        parent1,parent2 = random.sample(parents,2)
        
        if adaptive==True:
                        
            fp = np.max([parent1.fitness,parent2.fitness])
            crossover_rate = kwargs["k1"]*(kwargs["fmax"]-fp)/(kwargs["fmax"]-kwargs["fmean"]) if fp>=kwargs["fmean"] else kwargs["k1"] 
            
        if np.random.rand()<=crossover_rate: #if crossover rate is high, high change of crossover
            
            new,ID = mating(parent1,parent2,mode,ID)
            offspring+=new
            
        else:
            # we have to deepcopy parameters of parents otherwise continuous string objects are copied in memory
            parameters_parent1 = parent1.deepcopy_parameters()
            parameters_parent2 = parent2.deepcopy_parameters()
            # deepcopy chromosome object            
            offspring1 = deepcopy(parent1);offspring1.setParameters(parameters_parent1)
            offspring2 = deepcopy(parent2);offspring2.setParameters(parameters_parent2)
            # and now add to offspring
            offspring += [offspring1,offspring2]
    
        if multi_objective==True:
        
            cond = len(offspring) < len(parents)
        
        else:
            cond = len(offspring)<n_chromosomes-len(parents)
        
    return parents,offspring,ID
    
def mating(parent1,parent2,mode,ID):
        
    # binary encoding for IVS
    if mode == "binary":
        offspring,ID = mating_binary(parent1,parent2,ID)

    # list of lists encoding for MI
    if mode == "variable":
        offspring,ID = mating_variable(parent1,parent2,ID)

    # continuous encoding for PE
    if mode == "continuous":
        offspring,ID = mating_continuous(parent1,parent2,ID)
            
    return offspring,ID
            
def mating_binary(parent1,parent2,ID):

    l = len(parent1.parameters)
    beta= int(np.round(np.random.rand(1)*l))    
    offspring = [0,0]
    parents = [parent1,parent2]
   
    for i in range(2):
        
        offspring[i] = deepcopy(parents[i])
        offspring[i].setID(deepcopy(ID));ID=ID+1
        offspring[i].setFitness(np.nan)
 
    # copy parameters for offspring 1
    parameters_offspring = parent1.parameters.copy()
    index = parameters_offspring.index

    parameters_offspring.loc[index[beta:],"sample"] = parent2.parameters.loc[index[beta:],"sample"]
    offspring[0] = chromosome(parameters_offspring);offspring[0].setID(deepcopy(ID));ID=ID+1
    
    # copy parameters for offspring 2
    parameters_offspring = parent1.parameters.copy()
    parameters_offspring.loc[index[:-(l-beta)],"sample"] = parent2.parameters.loc[index[:-(l-beta)],"sample"]
    offspring[1] = chromosome(parameters_offspring);offspring[1].setID(deepcopy(ID));ID=ID+1
#        
    return offspring,ID

def mating_variable(parent1,parent2,ID):
    
    betamax  = 2
    alpha = np.random.randint(0,len(parent1.parameters))
    offspring = [0,0]
    parents = [parent1,parent2]
   
    for i in range(2):
        
        offspring[i] = deepcopy(parents[i])
        offspring[i].setID(deepcopy(ID));ID=ID+1
        offspring[i].setFitness(np.nan)
        
    sample1 = parent1.parameters["sample"]
    sample2 = parent2.parameters["sample"]
    
    # the next lines are actually chinese    
    "first offspring"
    index = sample1.index
    
    for i in range(alpha,len(index)):
        if (type(sample1.loc[index[i]])!=int) & (type(sample2.loc[index[i]])!=int):
            new = deepcopy(offspring[0].parameters.loc[index[i],"sample"])
            # crossover with beta weighting sample 1 and 2
            beta = np.random.uniform(0,betamax,len(sample1.loc[index[i]].returnString()))
            new.setString(sample1.loc[index[i]].returnString()-beta*(sample1.loc[index[i]].returnString()-sample2.loc[index[i]].returnString()))
            new.fixBoundary()
            offspring[0].parameters.loc[index[i],"sample"] = new
        else:
            offspring[0].parameters.loc[index[i],"sample"] = deepcopy(sample2.loc[index[i]])
 
    for i in range(alpha,len(index)):
        if (type(sample1.loc[index[i]])!=int) & (type(sample2.loc[index[i]])!=int):
            new = deepcopy(offspring[1].parameters.loc[index[i],"sample"])
            beta = np.random.uniform(0,betamax,len(sample1.loc[index[i]].returnString()))
            new.setString(sample2.loc[index[i]].returnString()+beta*(sample1.loc[index[i]].returnString()-sample2.loc[index[i]].returnString()))  
            new.fixBoundary()
            offspring[1].parameters.loc[index[i],"sample"] = new
    
        else:
            offspring[1].parameters.loc[index[i],"sample"] = deepcopy(sample1.loc[index[i]])

    return offspring,ID

def mating_continuous(parent1,parent2,ID):
    
    offspring = [0,0]
    parents = [parent1,parent2]
   
    for i in range(2):
        
        offspring[i] = deepcopy(parents[i])
        offspring[i].setID(deepcopy(ID));ID=ID+1
        offspring[i].setFitness(np.nan)
        
    sample1 = parent1.parameters["sample"]
    sample2 = parent2.parameters["sample"]

    beta = np.random.uniform(0,1.5)
    
    index = sample1.index
    
    for i in index:
        
        new1 = deepcopy(offspring[0].parameters.loc[i,"sample"])
        new2 = deepcopy(offspring[0].parameters.loc[i,"sample"])

        new1.setString(sample1.loc[i].returnString()-beta*(sample1.loc[i].returnString()-sample2.loc[i].returnString()))
        new1.fixBoundary()

        new2.setString(sample2.loc[i].returnString()+beta*(sample1.loc[i].returnString()-sample2.loc[i].returnString()))  
        new2.fixBoundary()
        
        offspring[0].parameters.loc[i,"sample"] = new1
        offspring[1].parameters.loc[i,"sample"] = new2

    return offspring,ID

def mutation(chromosomes,mutation_rate,mode,ID,adaptive=False,**kwargs):
    
    for i in range(len(chromosomes)):
        
        if chromosomes[i].protect == False:
            mutation_flag = chromosomes[i].mutation_operator(mutation_rate,mode,ID,adaptive=adaptive,**kwargs)
            if mutation_flag == True:
                ID = ID+ 1
        else:
            chromosomes[i].setProtect(False)
    return chromosomes,ID

def create_dir(res,L):
    
    for i in range(len(L)):
        if not os.path.exists(os.path.join(res,L[i])):
            os.makedirs(os.path.join(res,L[i]))      
    
def print_message(population,run,dt,multi_objective,res):
    
    if multi_objective==True:
        
        TSS = [i.performance["TSS"] for i in population.chromosomes if i.rank==1]
        print("Generation ("+str(run)+"): Max(TSS) = %.3f (average first front (+-std): %.3f (%.3f))"%(np.max(TSS),np.mean(TSS),np.std(TSS))+"\t (runtime: "+str(dt)+" seconds)")

    else:
        fitness = [i.fitness for i in population.chromosomes]
        print("Generation ("+str(run)+"): Max(fitness) = %.3f (average (+-std): %.3f (%.3f))"%(np.max(fitness),np.mean(fitness),np.std(fitness))+"\t (runtime: "+str(dt)+" seconds)")
      
    #write message to file
    if os.path.exists(os.path.join(res,"iterations.txt")):
        mode = "a"
    else:
        mode = "w"
    f = open(os.path.join(res,"iterations.txt"),mode)
    f.write("%i"%int(run)+","+str(population.best.fitness)+"\n")
    f.close()

    
