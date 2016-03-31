
Python scripts to develop and optimise a filter-based species distribution model

© 2016, Sacha Gobeyn (sacha.gobeyn@ugent.be or sachagobeyn@gmail.com). Licensed under CC BY 4.0 Creative Commons.

Details: Script can be run from a (I)python editor or command line
 
Publication: This code is coupled to a publication, it is advised to first read the publication.

Usage:

  * Model is run by running 'script.py' in Python OR in command line 'python script.py -flags'
   -flags: 

  * parameterfile
   -inpudata \t .csv file with inputdata

   -taxon                     name of taxon to develop and optimise model
 
   -variables                 .csv file with variables to consider in model

   -resmap                    name of map to print results
 
   -ga_settings               name of file containing options for genetic optimisation algorithm

   -md_settings               name of file containing options for model development


  * modelsettings
   -a1,a2,a3,a4               parameters of unimodal non-uniform habitat preference curves

  * gasettings (for detailed information, see Haupt and Haupt, 2004):
   -number_of_chromosomes     number of chromosomes to run GA

   -selection_rate            selection rate of population per generation
 
   -crossover_rate            crossover rate of population
 
   -mutation_rate             mutation rate of genes in chromosomes
 
   -objective_function        the objective function to minimize
 
   -maximum_runs              maximum number of generations to calculate
 
   -stop_criteria             number of generations for which no improvement is observed in the best value of the objective function.
 
   -ncpu                      number of cpu's to use on machine (-1 = all)
 
   -mode                      binary (a continues varient is implemented but not properly tested)

   -full_output               print all outputs of every chromosoom tested
 
   -nan_value                 nan-value for objective and evalutation criteria for infeasible solutions (has to be a number, int or float). Typically set way above the values of the -to minize!- objective function

Requirements:
  * Python 2.7.10
  
  * Numpy 1.9.2  
  
  * Pandas 0.16.2

  
  (tested on Python 2.7.10 with Anaconda 2.3.0 (64-bit)) 
  (Anaconda: free Python distribution package: https://www.continuum.io/downloads)
  

Support: Please feel free to open an issue in case of a problem.

References:
 Haupt,R.L., Haupt, S.E., 2004, Algorithms Practical Genetic Algorithms. John Wiley & Sons, Inc., Hoboken, New Jersey, 2nd edition, 251 pp.
