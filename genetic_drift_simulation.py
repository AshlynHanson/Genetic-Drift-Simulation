# ---------------------------------------------------------------------------------------------------------------
# Main purpose of the program:
# Molecular Evolution Simulation - Speciation Event (Genetic Drift)
#
# This program simulates the molecular evolution of a population during and after a speciation event. 
#
# The population will start with N haploid individuals of length, L. For the first couple of generations, 
# the population undergoes neutral evolution with the constant size N and a mutation rate of μ per position 
# per generation. At generation, g, the population will split into two equal populations of size, N/2 and evolution 
# will proceed for those two individual populations. This evolution process will continue until generation, G.
#
# Instructions:
# When run the program will prompt the user to enter the following information:
#
#   Population size (N)      The number of haploid individuals in the population (an even integer greater than 0)
#   Sequence length (L)      The number of base pairs in the sequence (an integer greater than 0)
#   Mutation rate (μ)        The rate at which to mutate the population (a decimal between 0 and 1)
#   Split generation (g)     The generation at which to split into two populations (an integer greater than 0)
#   Total Generations (G)    The total number of generations to run the simulation with (an integer greater than 0)
#
# The user will be prompted to choose whether or not to have every generation genotype displayed. If the user 
# chooses yes (Y), every genotype will be printed out for every generation. If the user chooses no (N), only the
# final generation genotype population will be displayed. This would be one population of size N or two populations
# of size N/2 depending on whether a split occurred at generation, g.

# Using the given information, the program will run a simulation of the speciation event, reproducing and mutating
# according to the mutation rate. The average genetic distance will be calculated for every generation and at the 
# end of the simulation will be displayed.
#
# Notes:
# If the inputted g is larger than the G inputted, no speciation event will occur. When the speciation event occurs 
# at g, the population is split at the g, meaning that every generation after that g value will be inside one of 
# two new populations. If the inputted population size, N is not an even number, the split populations will be 
# different sizes. The first population will have one more sequence than the second.
# 
# It is assumed that the population size, length of the sequence, and generations, g and G, will be greater 
# than or equal to 0. The mutation rate should be a decimal value between 0 and 1 inclusively.
#
# (C) 2021 Ashlyn Hanson, Yujun Zhi
# email: adh5584@truman.edu, yz4586@truman.edu
# ---------------------------------------------------------------------------------------------------------------

import random


def initialize_original_sequence(population_size_N, num_of_base_pairs_L):
    '''
    This function initializes the original generation sequence with the the 
    value 'A'.

    :param int population_size_N: the number of haploid individuals in the population
    :param int num_of_base_pairs_L: the length of base pairs for the sequence
    :return list initial_sequence: the original sequence of all adenine, 'A's
    '''
    initial_sequence = []
    for index in range(0, population_size_N):
        indiv_sequence = ""
        for nucleotide in range(0, num_of_base_pairs_L):
            indiv_sequence += "A"
        initial_sequence.append(indiv_sequence)
    return initial_sequence


def reproduce(mutation_rate, prev_generation):
    '''
    This function generates every new generation of the DNA sequences based on the previous generation.
    The mutate function is called within every generation to ensure that all new sequences are mutated
    according to the mutation rate. The mutation rate is expected to be a decimal value between 0 and 1

    :param float mutation_rate: the mutation rate per position per generation (between 0 and 1)
    :param list prev_generation: the previous generation of sequences to reproduce with mutation
    :return list new_generation: a list of sequences in the newly mutated population
    '''
    new_generation = []
    for index in range(0, len(prev_generation)):
		# mutate each sequence in the list of each individual sequence in the generation
        new_generation.append(mutate(mutation_rate, prev_generation[index]))
    return new_generation


def mutate(mutation_rate, generation_sequence):
    '''
    This function mutates the sequence according to the mutation rate. A random number is generated to determine
    which individual nucleotide in the sequence may or may not be mutated. Using the random library, a decimal number is generated between
    0 and 1. If the random number is within a specific interval inside the mutation rate, that nucleotide is mutated. To determine
    which nucleotide it becomes, a random choice is made from a list of C, T, and G. Because the initial sequence is assumed to be A,
    there are no mutations to Adenine or A. 

    :param float mutation_rate: the decimal mutation rate between 0 and 1 to determine how much the sequence will be mutated
    :param str generation_sequence: the previous generation sequence to mutate
    :return str new_generation_sequence: the new string generation sequence to be added to the new population
    '''
    new_generation_sequence = ""
    nucleotide_choices = ['C', 'T', 'G']
    for index in range(0, len(generation_sequence)):
        likelihood = random.random() # get a random floating point number between 0 and 1
        # if the random number is in the mutation rate percentage, mutate the nucleotides
        if likelihood <= mutation_rate:
            new_generation_sequence += random.choice(nucleotide_choices)
        # if the random number is outside the percentage interval, keep the nucleotide the same
        elif likelihood > mutation_rate:
            new_generation_sequence += generation_sequence[index]
    return new_generation_sequence


def genetic_distance(generation, population_size_N, num_of_base_pairs_L):
    '''
    This function calculates the average proportion of base pair differences between individuals within a 
    population. The total distance between each base pair is calculated for every haploid sequence in each
    individual population. These distances are then added together and the average is calculated by dividing by the 
    population size. Then, the average distance is calculated by dividing that by the length of the base pairs.

    For example, if a population consists of three individuals who differ from each other by 8, 10, and 14 
    base pairs in a sequence of length 100 bp, the average genetic distance within this population would be 
    (8 + 10 + 14) / 3 × (1/100), or 0.1067.

    :param int generation: the generation number from 1 to the population size (N)
    :param int population_size_N: The number of haploid individuals in the population
    :param int num_of_base_pairs_L: the length of the base pairs in the sequence
    :return float average_distance: The average genetic distance for this generation
    '''
    # calculate the distance between 2 sequences
    distances = []
    for index in range(0, len(generation)):
        for index2 in range(index + 1, len(generation)):
            distances.append(compare_sequences(generation[index], generation[index2]))
    
    # Add the sum of the distances together
    sum_of_distances = 0
    for value in distances:
        sum_of_distances += value

    average_distance = (sum_of_distances / population_size_N) / num_of_base_pairs_L

    return average_distance
	

def compare_sequences(sequence1, sequence2):
    '''
    This function compares two sequences and returns the number of different chars between the two.
    For this program, this function would compare two sequences in a generation of a population and
    returns the distance between the two. It is assumed that the two sequences are the same length.

    :param str sequence1: the first string of nucleotides to compare
    :param str sequence2: the second string of nucleotides to compare
    :return int num_differences: the number of differences between the two sequences, the distance
    '''
    num_differences = 0
    for index in range(0, len(sequence1)):
        if sequence1[index] != sequence2[index]:
            num_differences += 1
    return num_differences


def split_populations(generation_sequence, population_size_N):
    '''
    This function splits the population into two smaller sequences of size N/2.
    The first sequence will contain the first half of the previous generation's
    N values. The second will contain the second half of the previous generation's N values.
    It is assumed that the population size, N will be even. If it isn't the first sequence
    will contain one more haploid sequence than the second.

    :param list generation sequence: the previous generation sequence to split
    :param int population_size_N: the population size of the population generation
    '''
    population1= []
    population2 = []
    for index in range(0, len(generation_sequence)):
        if index < population_size_N / 2:
            population1.append(generation_sequence[index])
        else:
            population2.append(generation_sequence[index])
    return population1, population2


def print_output(avg_distance, generations):
    '''
    This function displays the average genetic distance per generation to the screen.
    The output would look similar to the example below:

    Generation      Avg Genetic Distance
    [  1            0.0                   ]
    [  2            6.3795                ]
    [  3            7.144500000000001     ]
    [  4            3.319                 ]
    [  4            3.22                  ]
    [  5            3.319                 ]
    [  5            3.22                  ]
    [  6            3.319                 ]
    [  6            3.22                  ]

    :param list avg_distance: a list of the average distances for every generation
    :param list generations: a list of all the number of generations run
    '''
    print('{:15s} {:20s}'.format('Generation', "Avg Genetic Distance"))
    for index in range(0, len(generations)):
        print('{:2s} {:12s} {:20s} {:2s}'.format("[ ", str(generations[index]), 
            str(avg_distance[index]), " ]"))


def display(sequence):
    '''
    This function prints the contents of an array, or an individual sequence, to the screen

    :param list sequence: the sequence to be displayed to the screen
    '''
    for index in range(0, len(sequence)):
        print("\t" + str(sequence[index]))


def print_instructions():
    '''
    This function prints the general program overview and instructions to the user
    '''
    print("\n===========================================================================================================")
    print('{:30s} {:50s}'.format(" ","Molecular Evolution Simulation - Speciation Event"))
    print("===========================================================================================================\n")

    print("PROGRAM OVERVIEW:\n")
    print("This program simulates the molecular evolution of a population during and after a speciation event.\n")

    print("The population will start with N haploid individuals or length, L. For the first couple of generations, the " +
        "population undergoes neutral evolution with the constant size N and a mutation rate of μ per position per generation. " +
        "At generation, g, the population will split into two equal populations of size, N/2 and evolution will proceed for " +
        "those two individual populations. This evolution process will continue until generation, G.\n")

    print("DNA STRUCTURE:\n") 
    print("The DNA sequence will contain only Adenine(A), Thymine(T), Cytosine(C), or Guanine(G). " +
        "Nucelotide ambiguity nucleotides will not be accepted for this program.\n")

    print("INPUT:\n")
    print("You will be asked to enter data for the following: \n")

    print('{:25s} {:50s}'.format('\tPopulation size (N)', 'The number of haploid individuals in the population (an even integer greater than 0)'))
    print('{:25s} {:50s}'.format('\tSequence length (L)', 'The number of base pairs in the sequence (an integer greater than 0)'))
    print('{:25s} {:50s}'.format('\tMutation rate (μ)', 'The rate at which to mutate the population (an decimal between 0 and 1)'))
    print('{:25s} {:50s}'.format('\tSplit generation (g)', 'The generation at which to split into two populations (an integer greater than 0)'))
    print('{:25s} {:50s}'.format('\tTotal Generations (G)', 'The total number of generations to run the simulation with (an integer greater than 0)\n'))

    print("INSTRUCTIONS:\n")
    print("Enter the necessary input information with the correct values. The simulation will then run with the given values and " +
        "output the average genetic distance for each population for each generation. NOTE: If the inputted g is greater than the inputted " +
        "G, the speciation event will not occur and there will not be a split.\n")
    print("Enter the requested data when prompted below to begin. \n")


def main():
    '''
    The main function prints the program overview and instructions to the user and then prompts the user
    for the necessary input to run the simulation. The inputted information is checked to ensure that all
    input is biologically possible. All sizes, lengths, and locations should be greater than 0. The user also
    has the option to choose whether to print every genotype for every generation. If they choose (Y) to print
    every generation, the output will be every generation genotype sequence displayed. If they choose (N) not
    to print every generation, only the final generation's population genotype will be displayed. At the end of
    the output, the average genetic distance for each generation will be displayed.
    '''
    print_instructions()
    # Get the user input
    # Check that the input matches the expected values and biologically plausible values (N >0, ect…)
    population_size_N = input("Please enter a population size (N): ")
    while int(population_size_N) <= 0:
        population_size_N = input("Please enter a valid population size greater than 0 (N): ")

    num_of_base_pairs_L = input("Please enter length of base pairs (L): ")
    while int(num_of_base_pairs_L) <= 0:
        num_of_base_pairs_L = input("Please enter a valid length of base pairs greater than 0 (L): ")

    mutation_rate = input("Please enter mutation rate (between 0 and 1): ")
    while float(mutation_rate) < 0 or float(mutation_rate) > 1:
        mutation_rate = input("Please enter a valid mutation rate (between 0 and 1): ")

    generation_g = input("Please enter generation g: ")
    while int(generation_g) <= 0:
        generation_g = input("Please enter a valid generation g that is greater than 0: ")

    total_generations_G = input("Please enter num of generations: ")
    while int(total_generations_G) <= 0:
        total_generations_G = input("Please enter num of generations that is less than generation g and greater than 0: ")

    list_all_generations = input("Would you like to see all of the generations' genotypes (Y or N): ")
    while list_all_generations.upper() != 'Y' and list_all_generations.upper() != 'N':
        list_all_generations = input("Please enter a valid input indicating yes or no (Y or N): ")

    # Convert the input from strings to the correct values
    population_size_N = int(population_size_N)
    num_of_base_pairs_L = int(num_of_base_pairs_L)
    mutation_rate = float(mutation_rate)
    generation_g = int(generation_g)
    total_generations_G = int(total_generations_G)

    # Initialize the initial generation with all ‘A’s
    initial_sequence = initialize_original_sequence(int(population_size_N), int(num_of_base_pairs_L))

    avg_genetic_distances = [] # The average genetic distances for each generation
    generations = [] # A list of all the generation numbers
    generation_number = 1 # Keeps track of the generation numbers

    if generation_g < total_generations_G:
        new_generation_sequence = initial_sequence
        for generation in range(0, generation_g):
            generations.append(generation_number)
            avg_genetic_distances.append(genetic_distance(new_generation_sequence, population_size_N, num_of_base_pairs_L))
            if list_all_generations.upper() == 'Y':
                print('-------- Generation: ' + str(generation_number) + ' --------\n')
                display(new_generation_sequence)
                print("Average genetic distance: " + str(genetic_distance(new_generation_sequence, population_size_N, num_of_base_pairs_L)) + "\n")
            new_generation_sequence = reproduce(mutation_rate, new_generation_sequence)
            generation_number += 1

        # Create 2 new positions with size N/2
        population1, population2 = split_populations(new_generation_sequence, population_size_N)
        for generation in range(generation_g, total_generations_G):
            generations.append(generation_number)
            avg_genetic_distances.append(genetic_distance(population1, population_size_N / 2, num_of_base_pairs_L))
            if list_all_generations.upper() == "Y":
                print('-------- Generation: ' + str(generation_number) + ' --------\n')
                print("Population 1: ")
                display(population1)
                print("Average genetic distance: " + str(genetic_distance(population1, population_size_N / 2, num_of_base_pairs_L)) + "\n")
            new_generation_sequence = reproduce(mutation_rate, population1)

            generations.append(generation_number)
            avg_genetic_distances.append(genetic_distance(population2, population_size_N / 2, num_of_base_pairs_L))
            if list_all_generations.upper() == "Y":
                print("Polulation 2: ")
                display(population2)
                print("Average genetic distance: " + str(genetic_distance(population2, population_size_N / 2, num_of_base_pairs_L)) +"\n")
            new_generation_sequence = reproduce(mutation_rate, population2)

            generation_number += 1
        if list_all_generations.upper() == "N":
                print('\n-------- Generation: ' + str(generation_number - 1) + ' --------\n')
                print("Population 1: ")
                display(population1)
                print("Average genetic distance: " + str(genetic_distance(population1, population_size_N / 2, num_of_base_pairs_L)) +"\n")
                print("\n")
                print("Polulation 2: ")
                display(population2)
                print("Average genetic distance: " + str(genetic_distance(population2, population_size_N / 2, num_of_base_pairs_L)) +"\n")

    # complete reproduction process from 0 to G
    else:
        generation_number = 1
        new_generation_sequence = initial_sequence
        for generation in range(0, total_generations_G):
            print('-------- Generation: ' + str(generation_number) + ' --------\n')
            generations.append(generation_number)
            if list_all_generations.upper() == "Y":
                display(new_generation_sequence)
            avg_genetic_distances.append(genetic_distance(new_generation_sequence, population_size_N, num_of_base_pairs_L))
            print("Average genetic distance: " + str(genetic_distance(new_generation_sequence, population_size_N, num_of_base_pairs_L)) +"\n")
            new_generation_sequence = reproduce(mutation_rate, new_generation_sequence)
            generation_number += 1
        if list_all_generations.upper() == "N":
            print('\n-------- Generation: ' + str(generation_number - 1) + ' --------\n')
            display(new_generation_sequence)
            print("Average genetic distance: " + str(genetic_distance(new_generation_sequence, population_size_N, num_of_base_pairs_L)) +"\n")

    print_output(avg_genetic_distances, generations)
    print("\n")


if __name__ == "__main__":
    '''
    Calls the main function to start the program
    '''
    main()
