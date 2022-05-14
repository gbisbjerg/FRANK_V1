from frank_helpers import *
from frank_exceptions import *
from ga_helpers import *
from gdc_helpers import *
from rgda_helpers import *
from initial_helpers import *
from alive_progress import alive_bar
import re , subprocess, os
        
def FRANK_header(settings):
    '''Looks cool
    https://ascii.co.uk/art/dinosaur
    '''
  
    print('''
                __
               / _)
        .-^^^-/ /
     __/       /    ''' + settings.name + '''
    <__.|_|-|_|
    ''')

def read_input(coordinatesFile, settingsFile):
    '''this function will take the user inputs and put it into whatever formats
       we actually want to use in the rest of the script. 
       Input: User inputs
       Outputs: ?? useful input, like initial geometry optimization card and
                settings saved for the single point calculations
       Uses helper functions'''
    settings = DefaultSettings()
    coordinates = coordinate_reader(coordinatesFile)
    geoInput, spInput, settings = settings_reader(settingsFile, settings)
    if settings.name == 'Unnamed Run':
        settings.name = coordinatesFile
    return coordinates, geoInput, spInput, settings

def initial_geo_op(coordinates, sampleGeoInput, settings):
    '''runs that initial geometry optimization because we don't trust the user lol
       Input: Gamess coordinates, geo input, and settings
       Output: nice structure, or we tell them the initial input was bad and 
               end the whole script'''
    log_string(settings, 'Initial Geometry Optimization')
    with alive_bar(1, bar=None, spinner="frank") as bar:
        bar.title('Initial Geometry Optimization')
        inputCard = settings.name + "/" + settings.name + 'initialGeoOp.inp'
        gamess_coords_to_gamess_input(settings, coordinates, sampleGeoInput, inputCard)

        outputFile = gamess_driver(settings, inputCard)
        try:
            error_check(outputFile)
            gamessCoords, bondMatrix = geo_op_output(outputFile, settings)

        except (SigsegvError, SigtermError, SCFConvergenceError, BugError):
            log_print(settings, 'A GAMESS error occurred. Restarting the failed calculation...')
            outputFile = gamess_driver(settings, inputCard)
            try:
                error_check(outputFile)
                gamessCoords, bondMatrix = geo_op_output(outputFile, settings)
            except SCFConvergenceError:
                log_print(settings, 'The SCF calculation failed to converge for the initial geometry optimization.')
                log_print(settings, 'Frank exiting...')
                quit()
                
            except:
                log_print(settings, 'Another GAMESS error occurred. This may be fixed by increasing memory in the geometry optimization input file.')
                log_print(settings, 'Frank exiting...')
                quit()
                
        except GamessMemoryError:
            log_print(settings, 'The GAMESS calculation failed due to insufficient memory. \nPlease increase the allocated memory in the geometry optimization input file and try again.')
            log_print(settings, 'Frank exiting...')
            quit()
            
        except GamessError:
            something_wrong_with_gamess(outputFile)
    
        except:
            log_print(settings, 'An unexpected error occured during the initial geometry optimization.')
            log_print(settings, 'Frank exiting...')
            quit()
        bar()
    return gamessCoords, bondMatrix

def population_generator(gamessCoords, bondMatrix, settings):
    '''creates a population using the nice geometry optimization structure as
       the default structure. By default structure, what I mean is we'll be using
       dihedral angles, and when we want cartesian coordinates, we can use the dihedral
       angles to change the original "default" structure. I feel like that will simplify 
       our code, but I'm down to change that if it doesn't lol 
       Input: geo op structure
       Output: population, a lot of dihedral angles
       Uses helper functions'''
    dihedralKey, everythingAfterKey = dihedral_key_generation(bondMatrix, gamessCoords)
    log_print(settings, "\nInitializing population...")
    initialPopulation = init_population(settings, dihedralKey, everythingAfterKey, gamessCoords)  
    return initialPopulation, dihedralKey, everythingAfterKey

def rgda(population, bondMatrix, settings):
    '''Resultant Gradient Descent Algorithm. scoot scooting the population to be
       feasible, based on our settings. hopefully there's a good library we can use!
       Input: the population, and settings that include constraints
       Output: new population
       Uses helper functions'''
    step_size = settings.RGDA_step_size
    iterations = settings.RGDA_max_iterations
    with alive_bar(len(population), bar= 'bubbles', spinner="fishes") as bar:
        for i, structure in enumerate(population):
            xyzCoords = [atom[2:] for atom in structure.coordinates]
            xyzCoords_corrected = RGDA_one_structure(xyzCoords, bondMatrix, step_size, iterations)
            population[i].coordinates = xyzCoords_to_gamess(xyzCoords_corrected, structure.coordinates)
            bar()
    return population

def gga(population, bondMatrix, dihedralKey, everythingAfterKey, sampleSPInput, sampleGeoInput, settings):
    '''The big bad boi, gradient-based genetic algorithm. A.k.a. just a loop that 
       calls the other parts of the algorithm
       Input: population, settings that include constraints
       Output: the global solution, hopefully'''
    log_print(settings, '\n\n   GGA Loop #' + str(settings.GGA_current_iteration))
    print('..-. .-. .- -. -.-')
    print('Running GDC...')
    population = gdc(population, bondMatrix, sampleGeoInput, settings)
    log_print(settings, 'GDC Complete.\n\nRunning GA...')
    population = ga(population, dihedralKey, everythingAfterKey, sampleSPInput, settings)
    log_print(settings, 'GA Complete.\n\nRunning RGDA...')
    population = rgda(population, bondMatrix, settings)
    log_print(settings, 'RGDA Complete.\n')
    settings.GGA_current_iteration += 1
    return population

def gdc(population, bondMatrix, sampleGeoInput, settings):
    '''Gradient Descent with the fitness function as Constraints. Basically 
       GAMESS geometry optimization
       Input: population, geometry optimization settings
       Output: new population
       Uses helper functions'''
    i = 0
    with alive_bar(len(population), bar= 'bubbles', spinner="fishes") as bar:
        while i < len(population):
            inputCard = settings.name + "/" + settings.name + str(population[i].id) + '.inp'
            gamess_coords_to_gamess_input(settings, population[i].coordinates, sampleGeoInput, inputCard)
            outputFile = gamess_driver(settings, inputCard)
            try:
                error_check(outputFile)
                gamessCoords, newBondMatrix = geo_op_output(outputFile, settings)
                energy = get_energy(outputFile)
                noError = True
                    
            except (SigsegvError, SigtermError, SCFConvergenceError, BugError):
                log_print(settings, 'A GAMESS error occurred. Restarting the failed calculation...')
                outputFile = gamess_driver(settings, inputCard)
                try:
                    error_check(outputFile)
                    gamessCoords, newBondMatrix = geo_op_output(outputFile, settings)
                    energy = get_energy(outputFile)
                    noError = True
                  
                except SCFConvergenceError:
                    log_print(settings, 'The SCF calculation failed to converge. Removing the following structure from the population:')
                    log_print(settings, 'Structure #' + population[i].id)
                    for j in range(len(population[i].coordinates)):
                        print(population[i].coordinates[j][0], population[i].coordinates[j][1], population[i].coordinates[j][2], population[i].coordinates[j][3], population[i].coordinates[j][4])
                    population.pop(i)
                    noError = False
                    
                except:
                    log_print(settings, 'Another GAMESS error occurred. If this is occuring often, please consider increasing the memory in the geometry optimization file.')
                    log_print(settings, 'Removing the following structure from the population:\nStructure #' + population[i].id)
                    for j in range(len(population[i].coordinates)):
                        print(population[i].coordinates[j][0], population[i].coordinates[j][1], population[i].coordinates[j][2], population[i].coordinates[j][3], population[i].coordinates[j][4])
                    population.pop(i)
                    noError = False
                    
            except GamessMemoryError:
                log_print(settings, 'The GAMESS calculation failed due to insufficient memory. \nPlease increase the allocated memory in the geometry optimization input file and try again.')
                log_print(settings, 'Frank exiting...')
                quit()
                
            except GamessError:
                something_wrong_with_gamess(outputFile)
        
            except:
                log_print(settings, 'An unexpected error occured during the GDC with GAMESS.')
                log_print(settings, 'Frank exiting...')
                quit()
            bar()
            if noError == True:
                if newBondMatrix == bondMatrix:
                    population[i].coordinates = gamessCoords
                    population[i].energy = energy
                    i += 1
                else:
                    population.pop(i)
		
    population = duplicate_check(population)        
    return population
        

def ga(population, dihedralKey, everythingAfterKey, sampleSPInput, settings):
    '''Genetic Algorithm. hopefully we can use a ga library. Our fitness function
       will be the single point calculation in GAMESS
       Input: population, settings for the GA
       Output: new population
       Uses helper functions'''
    #This loop iterates over the generations
    generations_best = [None] * settings.GA_generations
    
    for i in range(settings.GA_generations):
        with alive_bar(len(population), bar= 'bubbles', spinner="fishes") as bar:
		#Creates the next generation or initalizes the population
            if len(population) == 1:
                population = init_population(settings, dihedralKey, everythingAfterKey, population[0].coordinates) 
            elif len(population) > 1:
                population = next_generation(settings, population, dihedralKey, everythingAfterKey)
            else:
                raise Exception("GA - A starting structure must be supplied")
    
    		#Runs the next generation through the fitness function and determins their energy
            j = 0
            while j < len(population):
                inputCard = settings.name + "/" + settings.name + str(population[j].id) + '.inp'
                population[j].energy = fitness_function(settings, population[j].coordinates, inputCard, sampleSPInput)
                if population[j].energy == 0:
                    population.pop(j)
                else:
                    j += 1
                bar()

        #breaks out of the loop if the population has converged
        generations_best[i] = copy.deepcopy(best_structure(population))
        if is_converged(settings, generations_best, i, "GA"):
            log_print(settings, "GA - Population converged after {} generations".format(i))
            break
    return population

def is_the_gga_done_yet(population, settings, gga_best):
    '''This function will check if the user's settings for a finished solution are
       satisfied (ex. max # of loops for gga, or population "converged")
       Input: population, settings for the gga
       Output: True if not finished, False if finished. counter-intuitive but for the
               while loop
       Uses helper functions'''
    # needs convergence criteria and other stopping criteria (ex. max time)
    if is_converged(settings, gga_best, settings.GGA_current_iteration-2, "GGA"):
        log_print(settings, "GGA - Population converged after {} generations".format(settings.GGA_current_iteration - 1))
        return False
    elif settings.GGA_current_iteration - 1 < settings.GGA_max_iterations:
        return True
    else:
        return False

def prettify_result(ggaOutput):
    '''Gives a user friendly output
       Input: GGA result
       Output: Nice output
       Uses helper functions'''
    pass
