import random
import sys
sys.path.append('../') # what does this do?
from frank_helpers import *
from initial_helpers import dihedral_key_generation, FindDihedral # move to frank_helpers?
import numpy as np
import copy


#This function generates an intial population by mutating the intial structure
def init_population(settings, dihedralKey, everythingAfterKey, gamessCords):
    init_population = [None] * settings.GA_population_size
    initial_structure = StructureNode(gamessCords) 
    init_population[0] = initial_structure  #Store the inital structre

    for i in range(1, settings.GA_population_size):
        init_population[i] = dihedral_mutation(dihedralKey, everythingAfterKey, initial_structure)
    
    return init_population

def next_generation(settings, initial_generation, dihedralKey, everythingAfterKey):  
    next_generation = next_generation_tournament(settings, initial_generation, dihedralKey, everythingAfterKey)
    return next_generation

#Generates strucutres to account for the missing population by randomly selecting structures 
#from the initial generation and conducting a point mutation
def restore_population_size(settings, initial_generation, dihedralKey, everythingAfterKey):
    population_shortage = settings.GA_population_size - len(initial_generation)
    restore_population = [None] * population_shortage

    for i in range(0, population_shortage):
        random_parent = initial_generation[random.randrange(len(initial_generation))]
        restore_population[i] = dihedral_mutation(dihedralKey, everythingAfterKey, random_parent)
    return restore_population

#Rearranges array so that a structure cannont be paired with itself.
#NOTE: a maximum of two values for the same structure is expected
#Not con
def non_consecutive_ordering(ordering):
    if len(ordering) > 1:
        #If last two values are the same swap the last with the first
        if ordering[-1] == ordering[-2]:
            ordering[0], ordering[-1] = ordering[-1], ordering[0]

        #if a pair of values are the same then swap with the next group
        for i in range(0, len(ordering) - len(ordering) % 2, 2):
            try:
                if ordering[i] == ordering[i + 1]:
                    a = i + 1
                    b = i + 2
                    ordering[b], ordering[a] = ordering[a], ordering[b]
            except:
                pass
    return ordering

#generates an array of the pairs to be matched in the tournament,
#NOTE: values are spaced to avoid pairing with themselves
def tournament_pair_ordering(pop_to_fill):
    pair_ordering = list(range(0, pop_to_fill)) + list(range(0, pop_to_fill))
    random.shuffle(pair_ordering)
    pair_ordering = non_consecutive_ordering(pair_ordering)

    return pair_ordering

#Generate an array of parent pairs to be used for the tournamnet ordering, 
#Parents are selected be having a lower energy then the structure the tournament pairs them with
#NOTE: values are spaced to avoid pairing with themselves ex. [S1, S1] is invalid
def tournament_parent_ordering(settings, initial_generation, pop_size):
    pop_to_fill = settings.GA_population_size - pop_size
    pair_ordering = tournament_pair_ordering(pop_to_fill)

    parent_ordering = [None] * pop_to_fill
    for i in range(pop_to_fill):
        index = 2*i
        #structure1
        structure1_index = pair_ordering[index]
        structure1_energy = initial_generation[structure1_index].energy

        #structure2
        structure2_index = pair_ordering[index + 1]
        structure2_energy = initial_generation[structure2_index].energy
        parent_ordering[i] = tournament_win(structure1_index, structure1_energy, structure2_index, structure2_energy)

    parent_ordering = non_consecutive_ordering(parent_ordering)
    return parent_ordering

#Determins the winner of a tournament 
def tournament_win(index1, energy1, index2, energy2):
    if energy1 <= energy2:
        return index1
    else:
        return index2

#Rolls to see if a structure will undergo a point mutation
def mutation_roll(settings, dihedralKey, everythingAfterKey, structure):
    mutation_roll = random.uniform(0, 1)
    if mutation_roll <= settings.GA_point_mutation_rate:
        structure = dihedral_mutation(dihedralKey, everythingAfterKey, structure)
    return structure

#Selection option for generating the next generation
#An elite population is preserved
#All structures (elite included) are paired up randomly, the lower energy structure is selected 
#Pairs of winners are "breed" together. Mutation is applied evenly throughout all structures
def next_generation_tournament(settings, initial_generation, dihedralKey, everythingAfterKey):
    #Mutate based on any randomly selected strucutre to climb back to pop-setting size
    next_generation = restore_population_size(settings, initial_generation, dihedralKey, everythingAfterKey)

    #In place sorts the array based on energy (Ascending)
    initial_generation.sort(key=lambda x: x.energy, reverse=False)

    #Safeguard the best solutions
    elite_population = initial_generation[:settings.GA_elite_count]
    next_generation += elite_population

    #Randomize the selection order for the tournament
    parent_ordering = tournament_parent_ordering(settings, initial_generation, len(next_generation))
    ordering_tracker = 0

    #If there is an odd number of parents
    if len(parent_ordering) % 2 == 1:
        last_structure = initial_generation[parent_ordering[-1]]
        next_generation += [dihedral_mutation(dihedralKey, everythingAfterKey, last_structure)]

    #Loop to generate offspring through a tournament method
    while len(next_generation) < settings.GA_population_size:      
        parent1 = parent_ordering[ordering_tracker]
        parent2 = parent_ordering[ordering_tracker + 1]
        ordering_tracker += 2

        child1, child2 = crossover(dihedralKey, everythingAfterKey, initial_generation[parent1], initial_generation[parent2])
        child1 = mutation_roll(settings, dihedralKey, everythingAfterKey, child1)
        child2 = mutation_roll(settings, dihedralKey, everythingAfterKey, child2)

        next_generation += [child1]
        if len(next_generation) < settings.GA_population_size:
            next_generation += [child2]

    return next_generation

def next_generation_elitism_parents(settings, initial_generation, dihedralKey, everythingAfterKey):
    next_generation = [None] * settings.GA_population_size

    #In place sorts the array based on energy (Ascending)
    initial_generation.sort(key=lambda x: x.energy, reverse=False)

    #Seperates out the fittest individuals for use of crossover and mutations
    parents = initial_generation[:settings.GA_elite_count]

    i = len(parents)
    for i in range(i):
        next_generation[i] = parents[i]
   
    #Loop to generate offspring randomly based on the parenting pool
    while i < settings.GA_population_size:
        parent1 = random.randint(0, settings.GA_elite_count - 1)
        parent2 = random.randint(0, settings.GA_elite_count - 1)
        method_selection = random.randint(0, 100)
        if method_selection % 2 == 0:  ###Curently an even split between the two methods but this can be optimized
            next_generation[i] = dihedral_mutation(dihedralKey, everythingAfterKey, parents[parent1])
            i += 1
        else:
            next_generation[i], child2 = crossover(dihedralKey, everythingAfterKey, parents[parent1], parents[parent2])  #need to get parent 2
            if (settings.GA_population_size - i) > 1:
                next_generation[i + 1] = child2
            i += 2
    return next_generation

def fitness_function(settings, structure, fileName, sampleInput):
    gamess_coords_to_gamess_input(settings, structure, sampleInput, fileName)
    outputFile = gamess_driver(settings, fileName)
    try:
        error_check(outputFile)
        energy = get_energy(outputFile)
    
    except (SigsegvError, SigtermError, SCFConvergenceError, BugError):
        print('A GAMESS error occurred. Restarting the failed calculation...')
        outputFile = gamess_driver(settings, fileName)
        try:
            error_check(outputFile)
            energy = get_energy(outputFile)
          
        except SCFConvergenceError:
            print('The SCF calculation failed to converge. Removing the following structure from the population:')
            for i in range(len(structure)):
                print(structure[i][0], structure[i][1], structure[i][2], structure[i][3], structure[i][4])
            return 0
            
        except:
            print('Another GAMESS error occurred. If this is occuring often, please consider increasing the memory in the single point input file.')
            print('Removing the following structure from the population:')
            for i in range(len(structure)):
                print(structure[i][0], structure[i][1], structure[i][2], structure[i][3], structure[i][4])
            return 0
            
    except GamessMemoryError:
        print('The GAMESS calculation failed due to insufficient memory. \nPlease increase the allocated memory in the single point input file and try again.')
        print('Frank exiting...')
        quit()
        
    except GamessError:
        something_wrong_with_gamess(outputFile)

    except:
        print('An unexpected error occured during the GA with GAMESS.')
        print('Frank exiting...')
        quit()
    return energy

#This function conducts a two-point crossover and returns two new StructureNode
def crossover(dihedralKey, everythingAfterKey, parent1, parent2):
    crossover_pt1 = random.randrange(0, len(dihedralKey))
    crossover_pt2 = random.randrange(0, len(dihedralKey))

    if crossover_pt1 == crossover_pt2:
        locations_to_swap = [crossover_pt1]
    else:
        locations_to_swap = list(range(crossover_pt1, crossover_pt2 + 1, 1)) + list(range(crossover_pt1, crossover_pt2 - 1, -1))

    child1 = copy.deepcopy(parent1)
    child2 = copy.deepcopy(parent2)

    for loc in locations_to_swap:
        angle1 = FindDihedral(parent1.coordinates, [dihedralKey[loc]])[0]
        angle2 = FindDihedral(parent2.coordinates, [dihedralKey[loc]])[0]

        child1.coordinates = dihedral_angle_update(dihedralKey[loc], everythingAfterKey[loc], child1.coordinates, angle2) 
        child2.coordinates = dihedral_angle_update(dihedralKey[loc], everythingAfterKey[loc], child2.coordinates, angle1)

    return child1, child2

#Mutatues one of the dihedral angles and returns a new StructureNode
def dihedral_mutation(dihedralKey, everythingAfterKey, structure):
    #Select the angle to mutate and determin the new angle
    mutation_location = random.randrange(0, len(dihedralKey))
    mutation_angle = round(random.uniform(1, 359), 15)

    #Need to update coordiantes based on the rotation
    i = mutation_location
    updatedCoordinates = dihedral_angle_update(dihedralKey[i], everythingAfterKey[i], structure.coordinates, mutation_angle)
    newStructure = StructureNode(updatedCoordinates)
    return newStructure


#This function updates the dihedral angle 
def dihedral_angle_update(dihedral, everythingAfter, coords, newAngle): 
    
    angle = FindDihedral(coords, [dihedral])[0] - newAngle - 3.62611207820248 # fave magic number
    
    vectors = []
    for i in range(len(everythingAfter)):
        vectors.append([dihedral[1], everythingAfter[i]])
        
    axis = np.array([coords[dihedral[2] - 1][2] - coords[dihedral[1] - 1][2], coords[dihedral[2] - 1][3] - coords[dihedral[1] - 1][3], coords[dihedral[2] - 1][4] - coords[dihedral[1] - 1][4]]) / np.linalg.norm(np.array([coords[dihedral[2] - 1][2] - coords[dihedral[1] - 1][2], coords[dihedral[2] - 1][3] - coords[dihedral[1] - 1][3], coords[dihedral[2] - 1][4] - coords[dihedral[1] - 1][4]]))
    
    newVector = []
    for i in range(len(vectors)):
        vector = np.array([coords[vectors[i][1] - 1][2] - coords[vectors[i][0] - 1][2],\
                           coords[vectors[i][1] - 1][3] - coords[vectors[i][0] - 1][3],\
                           coords[vectors[i][1] - 1][4] - coords[vectors[i][0] - 1][4]])
    
        newVector.append((vector*np.cos(angle)) + (np.cross(axis, vector)*np.sin(angle)) + (axis*np.dot(axis, vector)*(1 - np.cos(angle))))
    
        coords[vectors[i][1] - 1][2] = coords[vectors[i][0] - 1][2] + newVector[i][0]
        coords[vectors[i][1] - 1][3] = coords[vectors[i][0] - 1][3] + newVector[i][1]
        coords[vectors[i][1] - 1][4] = coords[vectors[i][0] - 1][4] + newVector[i][2]        
        
    return coords