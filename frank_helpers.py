import re, os, sys, shutil, subprocess, platform, contextlib
from frank_exceptions import *

class DefaultSettings:
    log_file = ""

    '''DefaultSettings contains the run conditions for the GGA, 
    values are set to a default run setting which can be modifed based on the user input'''
    def __init__(self):
        self.GGA_max_iterations =  2
        self.GGA_convergence_threshold = 0.01 #inclusive
        self.energy_comparision = 5 # number of places after decimal
        self.GA_population_size = 5
        self.GA_convergence_threshold = 0.01 #inclusive
        self.GA_generations = 2
        self.GA_elite_count = 3
        self.GA_crossover_rate = 0.2
        self.GA_point_mutation_rate = 0.1
        self.RGDA_step_size = 0.001
        self.RGDA_max_iterations = 1000
        self.gamess_path = ""
        self.initial_bond_threshold_single = 0.5
        self.initial_bond_threshold_double = 1.5
        self.initial_bond_threshold_triple = 2.5
        self.name = 'Unnamed Run'
        self.cores = 1
        self.frank_path = ''
        self.GGA_current_iteration = 1 

    def __repr__(self):
        settings_string = \
        "Settings\n\
         GGA_max_iterations: %s\n\
         GGA_convergence_threshold: %s\n\
         energy_comparision: %s\n\
         GA_population_size: %s\n\
         GA_convergence_threshold: %s\n\
         GA_generations: %s\n\
         GA_elite_count: %s\n\
         GA_crossover_rate: %s\n\
         GA_point_mutation_rate: %s\n\
         RGDA_step_size: %s\n\
         RGDA_max_iterations: %s\n\
         gamess_path: %s\n\
         initial_bond_threshold_single: %s\n\
         initial_bond_threshold_double: %s\n\
         initial_bond_threshold_triple: %s\n\
         name: %s\n\
         cores: %s\n\
         frank_path: %s\n\
         GGA_current_iteration: %s\n\
         " % (self.GGA_max_iterations, self.GGA_convergence_threshold, self.energy_comparision, self.GA_population_size,  \
            self.GA_convergence_threshold, self.GA_generations, self.GA_elite_count, self.GA_crossover_rate, \
            self.GA_point_mutation_rate, self.RGDA_step_size, self.RGDA_max_iterations, self.gamess_path, \
            self.initial_bond_threshold_single, self.initial_bond_threshold_double, self.initial_bond_threshold_triple,\
            self.name, self.cores, self.frank_path, self.GGA_current_iteration)
        
        return settings_string

    def set_log_file(self):
        self.log_file = self.name + "/" + self.name + ".log"
        
class StructureNode:
    id_counter= 0
    '''StructureNode contains the energy and the XYZ coordinates for a 
       structure, and are the main base unit passed around'''
    def __init__(self, coordinates):
        self.id = StructureNode.id_counter
        self.energy = 0
        self.coordinates = coordinates
        StructureNode.id_counter += 1
        
    def __eq__(self, other):
        decimal_round = 8
        return ((type(other) == StructureNode)
          and round(self.energy, decimal_round) == round(other.energy, decimal_round)
        )
        
    def __le__(self, other):
        if self.energy <= other.energy:
            return True
        else:
            return False

    def __lt__(self, other):
        if self.energy < other.energy:
            return True
        else:
            return False

    def __repr__(self):
        c = self.coordinates

        string_c = ""
        for row in c:
            string_c += ("{: >4} {: >4} {: >15} {: >15}  {: >15} \n".format(*row))
        return "Structure %d\nEnergy: %d \n%s" % (self.id, self.energy, string_c) 

def gamess_coords_to_gamess_input(settings, structure, sampleInput, outInput):
    '''Takes the coordinates of the structure in GAMESS format, a sample
       GAMESS input card, and the file name for the new GAMESS input card 
       which the function writes'''
    inputCard = open(sampleInput, 'r')
    inputLines = inputCard.readlines()
    inputLines2 = inputCard.read()
    inputCard.close()
    outInput = open(outInput, 'w+')
    try:
        if re.findall('$DATA', inputLines2) != None:
            i = 0
            while True:
                '''Writes all the lines prior to $DATA from the sample input
                   card to the new input card'''
                if inputLines[i][:6] == ' $DATA':
                    break
                outInput.write(inputLines[i])
                i += 1
        else:
            raise ReallyBadInput
    except ReallyBadInput:
        print('Error: The sample input does not contain an $DATA group.')
        print('Exiting...')
        quit()
    outInput.write(' $DATA\n')
    outInput.write('Current Structure \n')
    outInput.write('C1 1\n')
    for i in range(len(structure)):
        outInput.write(' '.join(str(x) for x in structure[i]) + '\n')
    outInput.write(' $END')
    outInput.close()

    log_geo_op_input(settings, structure)


def clear_dat_files():
    files = os.listdir('.')
    fileCount = 0

    for file in files:
        end = file[-4:]
        if end == '.dat':
            os.remove(file)
            fileCount += 1

def gamess_driver( settings, inputCard ):
    my_os = platform.system()

    if my_os == "Windows":
        return gamess_driver_windows(settings, inputCard)
    elif my_os == "Darwin": #Apple
        return gamess_driver_mac(settings, inputCard)
    elif my_os == "Linus":
        pass
    else:
        raise Exception("GAMESS Driver - Unrecognized Operating System")

def gamess_driver_mac( settings, inputCard ):
    gamessPath = settings.gamess_path
    frankPath = settings.frank_path

    inputCard_noEnd = inputCard[:-4]
    cmd = gamessPath + " " + inputCard
    p = subprocess.Popen( cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8' )
    p.communicate(input= frankPath + inputCard_noEnd + ".log", timeout=None)
    clear_dat_files()
    return inputCard_noEnd + '.log'

def gamess_driver_windows(settings, inputCard, frankPath):
    gamessPath, numCores = settings.gamess_path, settings.cores
    frankPath = settings.frank_path

    '''Makes calls to GAMESS with the previously generated input cards'''
    inputCard_noEnd = inputCard[:-4]
    cmd = 'cd ' + gamessPath + '&rungms.bat ' + frankPath + inputCard + ' 2020.R2.pgiblas ' + str(numCores) + ' ' + frankPath + inputCard_noEnd + ".log"
    p = subprocess.call(cmd, shell=True)
    return inputCard_noEnd + '.log'

def geo_op_output(outputFile, settings):
    '''Reads the output of a geometry optimization and returns the gamess coordinates
       in an array and the bond order array'''

    bondThresholdSingle = settings.initial_bond_threshold_single
    bondThresholdDouble = settings.initial_bond_threshold_double
    bondThresholdTriple = settings.initial_bond_threshold_triple

    
    output = open(outputFile, 'r')
    rawOutput = output.read()
    output.close()
    gamessMatch = re.findall("EQUILIBRIUM GEOMETRY LOCATED.*?INTERNUCLEAR", rawOutput, re.DOTALL)
    gamessOutput = gamessMatch[0].split('\n')
    xyzCoords = []
    for i in range(len(gamessOutput)):
        if 3 < i < (len(gamessOutput) - 2):
            outputRow = gamessOutput[i].split(' ')
            noSpacesRow = [x for x in outputRow if x]
            xyzCoord = []
            for j in range(len(noSpacesRow)):
                if j == 0:
                    xyzCoord.append(noSpacesRow[j])
                elif j >= 1:
                    xyzCoord.append(float(noSpacesRow[j]))
            xyzCoords.append(xyzCoord)
    
    
    bondMatch = re.findall("ATOM PAIR DIST  ORDER.*?TOTAL", rawOutput, re.DOTALL)
    bondOutput = bondMatch[len(bondMatch) - 1].split('\n')
    bondOrder = [None] * len(xyzCoords)
    for i in range(len(bondOrder)):
        bondOrder[i] = [None] * len(xyzCoords)
    
    for i in range(len(bondOutput)):
        if 0 < i < (len(bondOutput) - 2):
            currentRow = bondOutput[i].split(' ')
            rowNumbers = [float(x) for x in currentRow if x]
            if len(rowNumbers) < 12: 
                if len(rowNumbers) == 4:
                    bondOrder[int(rowNumbers[1]) - 1][int(rowNumbers[0]) - 1] = rowNumbers[3]
                else:
                    bondOrder[int(rowNumbers[1]) - 1][int(rowNumbers[0]) - 1] = rowNumbers[3]
                    bondOrder[int(rowNumbers[5]) - 1][int(rowNumbers[4]) - 1] = rowNumbers[7]
            else:
                bondOrder[int(rowNumbers[1]) - 1][int(rowNumbers[0]) - 1] = rowNumbers[3]
                bondOrder[int(rowNumbers[5]) - 1][int(rowNumbers[4]) - 1] = rowNumbers[7]
                bondOrder[int(rowNumbers[9]) - 1][int(rowNumbers[8]) - 1] = rowNumbers[11]
                
             
    for i in range(len(bondOrder)):
        for j in range(len(bondOrder[i])):
            if j < i:
                if bondOrder[i][j] == None:
                    bondOrder[i][j] = 0
                elif bondOrder[i][j] > bondThresholdTriple:
                    bondOrder[i][j] = 3
                elif bondOrder[i][j] > bondThresholdDouble:
                    bondOrder[i][j] = 2
                elif bondOrder[i][j] > bondThresholdSingle:
                    bondOrder[i][j] = 1
                else:
                    bondOrder[i][j] = 0
                    
    log_geo_op_output(settings, xyzCoords)
    return xyzCoords, bondOrder

def get_energy(outputFile):
    '''Takes the output file of a successfulS GAMESS run, and either returns the
       total energy in Hartrees or raises an error if the SCF did not converge'''
    output = open(outputFile, 'r')
    rawOutput = output.read()
    output.close()
    energyMatch = re.search('TOTAL ENERGY = ....(\d+\.\d+)', rawOutput)
    energyString = energyMatch[0]
    energyArray = energyString.split(' ')
    energy = float(energyArray[len(energyArray) - 1])
    if energy == 0:
        raise SCFConvergenceError
    else:
        return energy

def print_full_atom(fullAtom):
	for atom in fullAtom:
		print('{0:<55} {1:<15} {2:<8}'.format(str(atom[0]),str(atom[1]),str(atom[2]) ))

#Retuns a boolean if convergence has occured
def is_converged(settings, generations_best, generation_number, loop_type):
    #Need multiple generations to determin convergence
    if generation_number < 1:
        return False

    if loop_type == "GA" and \
    abs(generations_best[generation_number - 1].energy - generations_best[generation_number].energy) <= settings.GA_convergence_threshold:
        return True
    elif loop_type == "GGA" and \
    abs(generations_best[generation_number - 1].energy - generations_best[generation_number].energy) <= settings.GGA_convergence_threshold:
        return True
    else:
        return False

#Returns the structure from the array with the lowest energy
def best_structure(structure_array):
    if len(structure_array) == 0:
        raise Exception("GA - No structures can be found in the current generation")

    structure_array.sort(key=lambda x: x.energy, reverse=False)
    return structure_array[0]

def newFolder( fileName ):
    '''Generates a new file with the provided name, NOTE any previous file with the name is removed'''
    try:
        os.makedirs(fileName)
        open(fileName + "/" + fileName + ".log", 'a').close()    
        print("Directory " + fileName +  " Created")
    except:
        shutil.rmtree(fileName)
        os.makedirs(fileName)
        open(fileName + "/" + fileName + ".log", 'a').close()        
        print( "Directory " + fileName +  " was replaced")
        
def print_coordiantes(coordinates):
    print("Starting Coordinates")
    
    for atom in coordinates:
        print("{: >4} {: >4} {: >15} {: >15}  {: >15}".format(*atom))
    print("\n")

def logSettings(settings, coordinates, sampleInputGeo, sampleInputSP):
    settings.set_log_file()
    log_file = settings.log_file

    with open(log_file, 'a') as log_:
        with contextlib.redirect_stdout(log_):
            print(settings, "\n")
            print_coordiantes(coordinates)

    with open(sampleInputGeo,'r') as InputGeo, open(log_file,'a') as log_:
        log_.write("Geometry Optimization Sample File\n")
        for line in InputGeo:
            log_.write(line)

    with open(sampleInputSP,'r') as InputSP, open(log_file,'a') as log_:
        log_.write("\nSingle Point Sample File\n")
        for line in InputSP:
            log_.write(line)

#Called in gamess_coords_to_gamess_input and logs coordinates run BEFORE geometry optimization
def log_geo_op_input(settings, xyzCoords):
    log_file = settings.log_file
    
    with open(log_file, 'a') as log_:
        log_.write("Structure SENT to GAMESS\n")
        for atom in xyzCoords:
            log_.write("{: >4} {: >4} {: >15} {: >15} {: >15}\n".format(*atom))
        log_.write("\n")

#Called in geo_op_output and logs coordinates AFTER the geometry optimization
def log_geo_op_output(settings, xyzCoords):
    log_file = settings.log_file
    
    with open(log_file, 'a') as log_:
        log_.write("Structure RETURNED from GAMESS\n")
        for atom in xyzCoords:
            log_.write("{: >4} {: >4} {: >15} {: >15} {: >15}\n".format(*atom))
        log_.write("\n")

def log_best_structure(settings, xyzCoords):
    log_file = settings.log_file
    
    with open(log_file, 'a') as log_:
        log_.write('''
                        __
                       / _)
                .-^^^-/ / ^
             __/       /    
            <__.|_|-|_| and the global solution is...
            \n''')
        for atom in xyzCoords:
            log_.write("{: >4} {: >4} {: >15} {: >15} {: >15}\n".format(*atom))
        log_.write("\n")

def log_print(settings, string):
    print(string)

    log_file = settings.log_file
    with open(log_file, 'a') as log_:
        log_.write(string+"\n")

def log_string(settings, string):    
    log_file = settings.log_file
    with open(log_file, 'a') as log_:
        log_.write(string+"\n")

def quit():
    '''Ends FRANK early if something is really wrong (bad sample input file)'''
    sys.exit()
	
def error_check(outputFile):
    '''Takes the output file of the GAMESS run, and raises errors if the
       calculation failed'''
    output = open(outputFile, 'r')
    rawOutput = output.read()
    output.close()
    if re.search('SIGSEGV', rawOutput) != None:
        print()
        raise SigsegvError
    elif re.search('SIGTERM', rawOutput) != None:
        raise SigtermError
    elif re.search('ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY', rawOutput) != None:
        raise GamessMemoryError
    elif re.search('The error is most likely to be in the application', rawOutput) != None:
        raise BugError
    elif (re.search('ILLEGAL', rawOutput) != None or 
          re.search('INCORRECT', rawOutput) != None or 
          re.search('EXECUTION OF GAMESS TERMINATED -ABNORMALLY-', rawOutput) != None or 
          re.search('DDI Process 0: error code 911', rawOutput) != None or 
          re.search('Fatal error detected', rawOutput) != None):
        raise GamessError
    else:
        return True
    
def something_wrong_with_gamess():
    '''Something unexpected went wrong with GAMESS, likely a mistake in the
       sample input file. This ends FRANK early'''
    output = open(outputFile, 'r')
    rawOutput = output.readlines()
    output.close()   
    print('Error: Something unexpected happened in the GAMESS run.\n')
    print('These are the last 15 lines of the GAMESS raw output:\n')
    i = len(rawOutput) - 16
    while i < len(rawOutput):
        print(rawOutput[i])
        i += 1
    print('\n\nFrank exiting...')
    quit()
