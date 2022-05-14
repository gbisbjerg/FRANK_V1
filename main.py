from frank_funcs import *
from frank_exceptions import *
from frank_helpers import *

def main():
    
    coordinateFile = sys.argv[1]
    settingsFile = sys.argv[2]

    coordinates, sampleInputGeo, sampleInputSP, settings = read_input(coordinateFile, settingsFile)
    
    FRANK_header(settings)
    
    newFolder(settings.name)
    logSettings(settings, coordinates, sampleInputGeo, sampleInputSP)
    gamessCoords, bondMatrix = initial_geo_op(coordinates, sampleInputGeo, settings)

    initialPop, dihedralKey, everythingAfterKey = population_generator(gamessCoords, bondMatrix, settings)
    
    print('\nRunning RGDA...')
    population = rgda(initialPop, bondMatrix, settings)
    print('RGDA Complete.\n')

    gga_best = [None] * settings.GGA_max_iterations
    while is_the_gga_done_yet(population, settings, gga_best):
        population = gga(population, bondMatrix, dihedralKey, everythingAfterKey, sampleInputSP, sampleInputGeo, settings)
        gga_best[settings.GGA_current_iteration - 2] = copy.deepcopy(best_structure(population))
            # gdc -> ga -> rgda

    print('''
                __
               / _)
        .-^^^-/ / ^
     __/       /    
    <__.|_|-|_| and the global solution is...
    ''')

    for i in range(len(population[0].coordinates)):
        print(population[0].coordinates[i][0], population[0].coordinates[i][1], population[0].coordinates[i][2], population[0].coordinates[i][3], population[0].coordinates[i][4])
    log_best_structure(settings, population[0].coordinates)

    #prettify_result(ggaOutput)


main()
