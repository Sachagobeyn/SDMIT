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

What can we do with SDMIT?
-------------------------
The species distribution model identification tool is a software tool based on machine learning that aims to train species distribution models with genetic algorithms. The implemented SDM is a habitat suitability model, aiming to define the relation 
of species response (here presence/absence) as a function of environmental gradients or feature. The genetic algorithm serves as an optimisation method performing:

1. Wrapped Feature selection: Searching for a set of input features which best describe observed presence/absence patterns (with parameters of species response curves have to be defined a priori).
                                                 
2. Embedded Feature Selection: Searching for a set of input features and parameter estimates which best describe observed presence/absence patterns.

3. Single or Multi-Objective Optimisation (SOO and MOO): Considering one (simple genetic algorithm) or more objective (non-dominated sorting genetic algorithm II) for the optimisation.
 
 Ok, which SDM are we optimising?
---------------------------------
The species distribution model currently implemented in the package is a simple SDM relating probability of species occurence to habitat suitability:

![PO](http://mathurl.com/ybwhlk2x.png)

with 

![HSI](http://mathurl.com/ybm8zpmt.png) 

and x, the input data for a number of features. The relation between SI and x is determined by a species response curve that is a gaussian-like topped of function ...

![fig0](docs/example-response_0.png)

... or a skewed topped of function

![fig1](docs/example-response_1.png)

In the encoding of the algorithm (see *settingsfile.txt*), this logistic increasing and decreasing functions before and after the optimal range is indicated with the boolean **logit**. Turn logit of and we get linear curves. In mathematical terms the logit function is described by:

![SI](http://mathurl.com/ycttr8xu.png)

with 

![p](http://mathurl.com/yc4m8bex.png)

Now the values of ![theta](http://mathurl.com/y74s5qu3.png) representing the optimal and optimal range (see two example figures above). Now the values of SI for different features have to be combined to one HSI value and this can be done by choosing different aggregation functions (see string **interference** in *settingsfile.txt*). For one you could take the minimum (**minimum**) ...

![min](http://mathurl.com/yavq2h89.png)

... or the mean (**mean**) ... 

![mean](http://mathurl.com/ybtoaaa9.png)

... or the geometrix mean (**squaredproduct**), which could be considered as most adequate from a theoretical point of view ...

![gm](http://mathurl.com/ycfp2r8e.png)

To finally determine species occurence one applies a threshold on the HSI (string or float **threshold** in *settingsfile.txt*):

![HSIthreshold](http://mathurl.com/yb4dhlkh.png)

` So basically, the implemented species distribution model is a model relating habitat suitability to species occurence by means of a number of species response curves and a HSI threshold. What do we remember for the application of the model in SDMIT?
 1. Set **logit** to **True** in the *settingsfile.txt* file if you want logistic increasing and decreasing functions describing the suboptimal conditions for SI. If **False** then linear functions will be set.
  2. Select an interference/aggregation (**interference** in *settingsfile.txt*) function to compute SI to HSI (**squaredproduct** or **minimum** or **mean**) in the *settingsfile.txt*
 3. Select a value for the **threshold** in the *settingsfile.txt*. One can also decide to maximise (**max**) the threshold based on the TSS (i.e. thresholds will vary over the models being optimised).
 `


Now how does this optimisation work?
-------------------------------------
As mentioned above, these are a number of mode of actions. First you can choose if you want to do 'wrapper' or 'embedded' feature selection. The difference is quite easy, in the first one relies on the parameter estimation for ![theta](http://mathurl.com/y74s5qu3.png) and does not ask the algorithm to change them to search for a good model. In case of the second, one asks the algorithm to also change the values of these parameters, so to find 'better' or more optimal solutions. In general, one would prefer the second approach, as the feature search can be influenced by the set parameter values. However, if one is very confident about the parameters (with help of expert knowledge), one can decide to use 'wrapper' feature selection. Now, the method implemented to facilitate both ways of model learning is based on genetic algorithms. We will not get into the details of this exaclty works, and go further the practical aspect:

` Choose in the *settings.txt* the **mode** (string)
1. **variable**: Embedded feature selection
2. **binary**: Wrapper feature selection
That was easy, no?
`



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
Deb, K., Pratap, A., Agarwal, S., Meyarivan, T., 2002. A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE Trans. Evol. Comput. 6, 182–197.

Gibbs, M.S., Dandy, G.C., Maier, H.R., 2008. A genetic algorithm calibration method based on convergence due to genetic drift. Inf. Sci. (Ny). 178, 2857–2869.

Haupt, R.L., Haupt, S.E., 2004. Algorithms Practical Genetic Algorithms, 2nd ed. John Wiley & Sons, Inc., Hoboken.

Mouton, A.M., De Baets, B., Goethals, P.L.M., 2010. Ecological relevance of performance criteria for species distribution models. Ecol. Modell. 221,1995–2002.


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


Licensed under [ CC BY 4.0 Creative Commons](https://creativecommons.org/licenses/by/4.0/ "CC")
