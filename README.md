# SDMIT
Python scripts to develop and perform model identification

© 2017, Sacha Gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com). 

Details: Script can be run from a (I)python editor or command line. Designed for batch simulations
 
Publication: This code is coupled to a publication, it is advised to first read the publication.

Usage: see SCRIPTS/script.py

Requirements (tested with):
  * Python 3.5.4
  
  * Numpy 1.13.1  
  
  * Pandas 0.20.3
  
  * Sklearn 0.19.0
  
  (tested with Python 3.5.4 with [Anaconda 3 (64-bit)](https://www.continuum.io/downloads "Anaconda"))

Support: Please feel free to open an issue in case of a problem.

Section 1: What can we do with SDMIT?
-------------------------
The species distribution model identification tool is a software tool based on machine learning that aims to train species distribution models with genetic algorithms. The implemented SDM is a habitat suitability model, aiming to define the relation 
of species response (here presence/absence) as a function of environmental gradients or feature. The genetic algorithm serves as an optimisation method performing:

1. Wrapped Feature selection: Searching for a set of input features which best describe observed presence/absence patterns (with parameters of species response curves have to be defined a priori).
                                                 
2. Embedded Feature Selection: Searching for a set of input features and parameter estimates which best describe observed presence/absence patterns.

3. Single or Multi-Objective Optimisation (SOO and MOO): Considering one (simple genetic algorithm) or more objective (non-dominated sorting genetic algorithm II) for the optimisation.
 
Section 2: Ok, which SDM are we optimising?
---------------------------------
The species distribution model currently implemented in the package is a simple SDM relating probability of species occurence to habitat suitability:

![PO](http://mathurl.com/ybwhlk2x.png)

with 

![HSI](http://mathurl.com/ybm8zpmt.png) 

and x, the input data for a number of features. The relation between SI and x is determined by a species response curve that is a gaussian-like topped of function ...

![fig0](docs/example-response_0.png)

... or a skewed topped of function

![fig1](docs/example-response_1.png)

In the encoding of the algorithm (see *settings.txt*), this logistic increasing and decreasing functions before and after the optimal range is indicated with the boolean **logit**. Turn logit of and we get linear curves. In mathematical terms the logit function is described by:

![SI](http://mathurl.com/ycttr8xu.png)

with 

![p](http://mathurl.com/yc4m8bex.png)

Now the values of ![theta](http://mathurl.com/y74s5qu3.png) representing the optimal and optimal range (see two example figures above). Now the values of SI for different features have to be combined to one HSI value and this can be done by choosing different aggregation functions (see string **interference** in *settings.txt*). For one you could take the minimum (**minimum**) ...

![min](http://mathurl.com/yavq2h89.png)

... or the mean (**mean**) ... 

![mean](http://mathurl.com/ybtoaaa9.png)

... or the geometrix mean (**squaredproduct**), which could be considered as most adequate from a theoretical point of view ...

![gm](http://mathurl.com/ycfp2r8e.png)

To finally determine species occurence one applies a threshold on the HSI (string or float **threshold** in *settings.txt*):

![HSIthreshold](http://mathurl.com/yb4dhlkh.png)


>So basically, the implemented species distribution model is a model relating habitat suitability to species occurence by means of a number of species response curves and a HSI threshold. What do we remember for the application of the model in SDMIT?
 >1. Set **logit** to **True** in the *settings.txt* file if you want logistic increasing and decreasing functions describing the suboptimal conditions for SI. If **False** then linear functions will be set.
>2. Select an interference/aggregation (**interference** in *settings.txt*) function to compute SI to HSI (**squaredproduct** or **minimum** or **mean**) in the *settings.txt*
>3. Select a value for the **threshold** in the *settingsfile.txt*. One can also decide to maximise (**max**) the threshold based on the TSS (i.e. thresholds will vary over the models being optimised).

Section 3: Now what can we optimise and how?
---------------------------------
As mentioned above, these are a number of mode of actions. First you can choose if you want to do 'wrapper' or 'embedded' feature selection. The difference is quite easy, in the first one relies on the parameter estimation for ![theta](http://mathurl.com/y74s5qu3.png) and does not ask the algorithm to change them to search for a good model. In case of the second, one asks the algorithm to also change the values of these parameters, so to find 'better' or more optimal solutions. In general, one would prefer the second approach, as the feature search can be influenced by the set parameter values. However, if one is very confident about the parameters (with help of expert knowledge), one can decide to use 'wrapper' feature selection. Now, the method implemented to facilitate both ways of model learning is based on genetic algorithms. We will not get into the details of this exaclty works, and go further the practical aspect:


>Choose in the *settings.txt* the **mode** (string)
>1. **variable**: Embedded feature selection
>2. **binary**: Wrapper feature selection
>That was easy, no?

Now we have defined what we want to identify, however, we also need to define what we want to optimise. Typically, this is a objective function based on a accuracy measure. Many measures are available (Mouton *et al.*, 2010), however, the true skill statistic is an often used statistic:

![TSS](http://mathurl.com/yb6zzuwk.png)

with (Sn = correct estimation of species presence, Sp = correct estimation of species absence):

![SNSP](http://mathurl.com/y9s9xdrh.png)

with: 

![confusionmatrix](http://mathurl.com/ybpbre3c.png)

So we can decide to maximise on TSS (single objective optimisation) or we could also optimise on both. Now, how would the latter work? It would mean we would have to optimise two criteria; the correct estimation of species presence, on the one hand, and species absence, on the other. However, the objective of optimising the two measures can be conflicting: SDMs that estimate species presence prefectly fail to estimate absence well and vice versa. Now, to identify this trade-off between objectives, specific algorithms are developed some time ago. Here, we use the non-dominated sort genetic algorithm II (NSGA-II) to identify the trade-off between Sn and Sp. Also here, the indication whether we want to choose for SOO or MOO are in the *setting.txt* file.


>1. If you want to do SOO, go to the *settings.txt* file and set **multi_objective** to False. Otherwise (MOO) set it to True.
>2. Choose you objective function **objective_function** in the *settings.txt* file, for SOO, this can be the TSS, Kappa, AUC, AIC, BIC or CCI (see Mouton *et al.*, 2010 for formulation). For MOO, one has to define multiple objectives and delimit them with a comma (e.g. 'Sn,Sp'). The current impementation has been tested for two objectives, but in practice should also work for three to four objectives.

Finally, there are a number of settings we need to define before we start the algorithm. Now, we could just use standard values for these 'hyper parameters', however, depending on how much time you have, you might want to tune these a bit, since one could end up with non-satisfying results when stopping the search algorithm to soon. Basically, we need to calculate a few things:

1. The number of algorithm iteration cycles.
2. The population size PS, the number of solutions are evaluated in one algorithm iteration cycle.
3. The mutation rate pm, determing the rate of random perturbations in the solutions (e.g. change a parameter value of  species response curve to a new random value).
4. The selection rate, determining how much solutions of an algorithm iteration cycle are copy-pasted to the next cyle, which we typically set to 0.5.
5. The crossover rate pc, determining how the solutions are 'combined' to new solutions (compare it with the concept of inheritance  from parents to offspring) which we can set to 1.

Now about the last two we don't really need to worry, however the first three are quite important since they interact with each other. Luckily, a smart guy, Matthew Gibbs came up with some guidelines to determine these three hyper parameters and it turns out that they work quite well. During my research, I did a number of tests and found them to work well. The guidelines:

1. Determine how much time you want the algorithm to run. You can calculate this with how long a fitness evaluation (a single SDM run) takes (for instance 0.1 second) and the time you have available (30 minutes = 1800 seconds), i.e. determine function evaluations by dividing the computer time available by the average time to compute the fitness function. Thus function evaluation is then - in this example - equal to 18000.

2. Solve next equation to find PS with M = 3 and l equal to 4 multiplied by number of variables (*for embedded feature selection!!!*) OR equal to the number of variables (*for wrapper feature selection*) (for l is 10, approximately a value of 50 is found for PS):

![FE](http://mathurl.com/ycrp8uje.png)

3. The mutation rate pm is equal to 5 divided by PS (in this example 0.1).

4. The (minimum) number of algorithm iteration cycles would then be equal to FE/PS (approximately 360 for N = 50 and FE = 18000).

> Compute an approximate value for the hyper parameters by using the [hyperparameters.py](scripts/hyperparameters.py) function. Set selection rate equal to 0.5, crossover rate equal to 1 in the *settings.txt*.

Section 4: Now let's try to get the code running
-------------------------------------
The Python code is developed in a way it can run as an executable by typing in command window 'python ./scripts/script.py parameterfile.txt'. As a consequence, the user interface is run by adapting text or csv file. At first hand this might seem like a bit of a struggle, however, when trying to move towards machine learning and uncertainty analysis on high performance platforms, this is a great advantage.

Step 1: Prepare input files 

1.[input data file](inputdata_Baetidae.csv): a comma delimited file with columns ID, X, Y, date, taxon, abundance, variable and value. A coupled taxon-variable list format is used. X and Y can be the coordinates of the center of a pixel or the exact coordinates of the location. Make sure that the name of taxon and variable are the same in the files [considered variables](considered_variables.csv), [model parameters](parameters_Baetidae.csv)

2.[considered variables file](considered_variables.csv): a comma delimited file with columns variable, consider, name_sim and unit. With the column consider, one can switch (value 1 or 0) variables on/off that we want to consider in the optimisation. The column 'name_sim' should carry the same name for the variables considere in [input data file](inputdata_Baetidae.csv) and [model parameters](parameters_Baetidae.csv) and should contain normal characters (a-z, 1-9, no ',''s).

3.[model parameters file](parameters_Baetidae.csv): a comma delimited file with columns taxon, value (leave blank), type (always continuous), b1 to b4, low, high and a1 to a1. The a1 and a4 are the values for ![theta](http://mathurl.com/y74s5qu3.png) and are an initial estimate by the user. The values for b1 to b4 are boundary values, with b1 and b2 bounding a1 and a2, and b3 and b4 bounding a3 and a4. Make sure the names of the variables match the column variable in [input data file](inputdata_Baetidae.csv) and name_sim in [considered variables file](considered_variables.csv)

4.[settings file](settings.txt): a file delimited by tabs indicating the settings for the model to run. A list of important settings are already reported above and all are summarized below:


| tag        | type | default   | explanation / notes |
| ---|:---:|:---:|
| nan_value      | 100000000000000 | nan value used in computation objective function|
| multi_objective     | False      | SOO (False) or MOO (True) |
| objective_function | TSS     |    objective function used to optimise. If multiple objectives are used, then one has to delimit the objectives with a comma (e.g. Sn,Sp) |
| ncpu      | 1 | number of cpu's used to run algorithm, for all use '-1' 
| number_of_chromosomes     | -      |   population size, see above (minimum eight) |
| maximum_runs | -      |    number of iterations |
| stop_criterion      | - | number of iterations for which not improvement in the objective function has been observed |
| mutation_rate     | -      | see above, between 0. and 1.   |
| crossover_rate | 1.0     |    see above, between 0. and 1. |
| selection_rate      | 0.5 | see above, between 0. and 1. |
| duplicates     | True     |   only tested for True |
| mode | variable     |    wrapper feature (binary) or embedded feature (variable) selection |
| logit      | True | logistic increasing and decreasing function (True) or linear (False) |
| interference    | squaredproduct      |  interference between SI values from response curves (see Section 1)  |

5.[code parameter file](parameterfile.txt): a file delimited by tabs linking files together so [code](scripts/script.py) can run in command line. 

References:
-----------
Deb, K., Pratap, A., Agarwal, S., Meyarivan, T., 2002. A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE Trans. Evol. Comput. 6, 182–197.

Gibbs, M.S., Dandy, G.C., Maier, H.R., 2008. A genetic algorithm calibration method based on convergence due to genetic drift. Inf. Sci. (Ny). 178, 2857–2869.

Haupt, R.L., Haupt, S.E., 2004. Algorithms Practical Genetic Algorithms, 2nd ed. John Wiley & Sons, Inc., Hoboken.

Mouton, A.M., De Baets, B., Goethals, P.L.M., 2010. Ecological relevance of performance criteria for species distribution models. Ecol. Modell. 221,1995–2002.

<!---
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
-->


Licensed under [ CC BY 4.0 Creative Commons](https://creativecommons.org/licenses/by/4.0/ "CC")
