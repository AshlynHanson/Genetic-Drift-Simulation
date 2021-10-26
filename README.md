# Genetic Drift Simulation

## Main purpose of the program:
Molecular Evolution Simulation - Speciation Event (Genetic Drift)

 This program simulates the molecular evolution of a population during and after a speciation event. 

The population will start with N haploid individuals of length, L. For the first couple of generations, 
the population undergoes neutral evolution with the constant size N and a mutation rate of μ per position 
per generation. At generation, g, the population will split into two equal populations of size, N/2 and evolution 
will proceed for those two individual populations. This evolution process will continue until generation, G.

## Instructions:
When run the program will prompt the user to enter the following information:

   Population size (N)      The number of haploid individuals in the population (an even integer greater than 0)
   Sequence length (L)      The number of base pairs in the sequence (an integer greater than 0)
   Mutation rate (μ)        The rate at which to mutate the population (a decimal between 0 and 1)
   Split generation (g)     The generation at which to split into two populations (an integer greater than 0)
   Total Generations (G)    The total number of generations to run the simulation with (an integer greater than 0)

 The user will be prompted to choose whether or not to have every generation genotype displayed. If the user 
 chooses yes (Y), every genotype will be printed out for every generation. If the user chooses no (N), only the
 final generation genotype population will be displayed. This would be one population of size N or two populations
 of size N/2 depending on whether a split occurred at generation, g.

 Using the given information, the program will run a simulation of the speciation event, reproducing and mutating
 according to the mutation rate. The average genetic distance will be calculated for every generation and at the 
 end of the simulation will be displayed.

## Notes:
 If the inputted g is larger than the G inputted, no speciation event will occur. When the speciation event occurs 
 at g, the population is split at the g, meaning that every generation after that g value will be inside one of 
 two new populations. If the inputted population size, N is not an even number, the split populations will be 
 different sizes. The first population will have one more sequence than the second.
 
 It is assumed that the population size, length of the sequence, and generations, g and G, will be greater 
 than or equal to 0. The mutation rate should be a decimal value between 0 and 1 inclusively.

 (C) 2021 Ashlyn Hanson, Yujun Zhi
 email: ashlyndhanson@gmail.com
