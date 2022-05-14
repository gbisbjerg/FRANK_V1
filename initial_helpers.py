import numpy as np
import math 
from copy import deepcopy
import queue
import configparser

import sys
sys.path.append('../')
from frank_helpers import *


#Class structure for network when recursivly tracking through bonds
class Node:
	def __init__(self, atom):
		self.atom = atom
		self.neighbors = []

	def add_neighbor(self, neighbor):
		self.neighbors.append(neighbor)

	def rm_neighbor(self, neighbor):
		self.neighbors.remove(neighbor)

	def display(self):
		print("Atom: " + str(self.atom))
		print("Neighbors: " + str(self.neighbors)[1:-1])

atomicNumbers = {'H':1, 'HE':2, 'LI':3, 'BE':4, 'B':5, 'C':6, 'N':7, 'O':8, 
                     'F':9, 'NE':10, 'NA':11, 'MG':12, 'AL':13, 'SI':14, 'P':15, 
                     'S':16, 'CL':17, 'AR':18, 'K':19, 'CA':20, 'SC':21, 'TI':22, 
                     'V':23, 'CR':24, 'MN':25, 'FE':26, 'CO':27, 'NI':28, 'CU':29,
                     'ZN':30, 'GA':31, 'GE':32, 'AS':33, 'SE':34, 'BR':35, 'KR':36,
                     'RB':37, 'SR':38, 'Y':39, 'ZR':40, 'NB':41, 'MO':42, 'TC':43,
                     'RU':44, 'RH':45, 'PD':46, 'AG':47, 'CD':48, 'IN':49, 'SN':50,
                     'SB':51, 'TE':52, 'I':53, 'XE':54, 'CS':55, 'BA':56, 'LA':57,
                     'CE':58, 'PR':59, 'ND':60, 'PM':61, 'SM':62, 'EU':63, 'GD':64,
                     'TB':65, 'DY':66, 'HO':67, 'ER':68, 'TM':69, 'YB':70, 'LU':71,
                     'HF':72, 'TA':73, 'W':74, 'RE':75, 'OS':76, 'IR':77, 'PT': 78,
                     'AU':79, 'HG': 80, 'TL':81, 'PB':82, 'BI':83, 'PO':84, 'AT':85,
                     'RN':86, 'FR':87, 'RA':88, 'AC':89, 'TH':90, 'PA':91, 'U':92,
                     'NP':93, 'PU':94, 'AM':95, 'CM':96, 'BK':97, 'CF':98, 'ES':99,
                     'FM':100, 'MD':101, 'NO':102, 'LR':103, 'RF':104, 'DB':105, 'SG':106,
                     'BH':107, 'HS':108, 'MT':109, 'DS':110, 'RG':111, 'CN':112, 'NH':113,
                     'FL':114, 'MC':115, 'LV':116, 'TS':117, 'OS':118}

def coordinate_reader(coordinateFile):
    file = open(coordinateFile, 'r')
    coordinates = []
    for line in file:
        lineWords = len(line.split())
        if lineWords > 1:
            split_strings = line.split()
            atom = split_strings[0]
            split_strings.insert(1, str(atomicNumbers.get(atom)))
            coordinate = [float(split_strings[i]) if i > 0 else split_strings[i] for i in range(len(split_strings))]
            coordinates.append(coordinate)
    file.close()        
    return coordinates

def settings_reader(settingsFile, settings):
    parser = configparser.ConfigParser()
    parser.read(settingsFile)
    if parser.has_section('Run'):
        if parser.has_option('Run', 'name'):
            settings.name = parser.get('Run', 'name')
        if parser.has_option('Run', 'cores'):
            settings.cores = int(parser.get('Run', 'cores'))
        if parser.has_option('Run', 'geo_input'):
            geoInput = parser.get('Run', 'geo_input')
        else:
            print('Error: The [Run] settings must contain the name of the sample Geometry Optimization input file (geo_input)')
            quit()
        if parser.has_option('Run', 'sp_input'):
            spInput = parser.get('Run', 'sp_input')
        else:
            print('Error: The [Run] settings must contain the name of the sample Single Point input file (sp_input)')
            quit()
    else:
        print('Error: The settings file must contain [Run] settings including:\n Name of the sample Geometry Optimization input file (geo_input)\n Name of the sample Single Point input file (sp_input) ')
        quit()
    
    if parser.has_section('Paths'):
        if parser.has_option('Paths', 'frank'):
            settings.frank_path = parser.get('Paths', 'frank')
        else: # Change to call to see the current path??
            print('Error: The [Paths] setting must include the Frank path (frank)')
            quit()
        if parser.has_option('Paths', 'gamess'):
            settings.gamess_path = parser.get('Paths', 'gamess')
        else:
            print('Error: The [Paths] settings must include the Gamess path (gamess)')
            quit()
    else:
        print('Error: The settings file must contain [Paths] settings including:\n Frank Path (frank)\n Gamess path (gamess) ')
        quit()
    
    if parser.has_section('GGA'):
        if parser.has_option('GGA', 'max_iterations'):
            settings.GGA_max_iterations = int(parser.get('GGA', 'max_iterations'))
    
    if parser.has_section('GA'):
        if parser.has_option('GA', 'pop_size'):
            settings.GA_population_size = int(parser.get('GA', 'pop_size'))
        if parser.has_option('GA', 'generations'):
            settings.GA_generations = int(parser.get('GA', 'generations'))
        if parser.has_option('GA', 'elite_count'):
            settings.GA_elite_count = int(parser.get('GA', 'elite_count'))
        if parser.has_option('GA', 'crossover'):
            settings.GA_crossover_rate = float(parser.get('GA', 'crossover'))
        if parser.has_option('GA', 'mutation'):
            settings.GA_point_mutation_rate = float(parser.get('GA', 'mutation'))
    
    if parser.has_section('RGDA'):
        if parser.has_option('RGDA', 'max_iterations'):
            settings.RGDA_max_iterations = int(parser.get('RGDA', 'max_iterations'))
        if parser.has_option('RGDA', 'step_size'):
            settings.RGDA_step_size = float(parser.get('RGDA', 'step_size'))
            
    return geoInput, spInput, settings
		
# #Identify if the provided atom in part of a cycle, and returns a cycle
###Note the search is terminated when one ring is found, ARD24091 only results in one ring
def atom_cycle_check(atom_number, updated_coords, original, path = [], edges = []):
	if atom_number in path and atom_number == original:
		return path
	for edge in updated_coords[atom_number-1][1]:
		if (atom_number < edge):
			current_edge = (atom_number, edge)
		else:
			current_edge = (edge, atom_number)

		if current_edge not in edges and current_edge[0] != current_edge [1]:
			updated_edge = deepcopy(edges)
			updated_edge.append(current_edge)
			updated_path = deepcopy(path)
			updated_path.append(atom_number)
			cycle_found = atom_cycle_check(edge, updated_coords, original, updated_path, updated_edge)
			if cycle_found != -1:
				return cycle_found
	return -1

#This function generates all the possible cycles found in the structure
def cycle_check(updated_coords):
	cycles = []
	for atom in range(1,len(updated_coords) + 1):
		cycle = atom_cycle_check(atom, updated_coords, atom, [], [])
		if (cycle != -1):
			cycles.append(cycle)

	unique_cycles = []
	unique_sorted_cycle = []
	for c in cycles:
		temp_c = deepcopy(c)
		temp_c.sort()
		if temp_c not in unique_sorted_cycle:
			unique_cycles.append(c)
			unique_sorted_cycle.append(temp_c)
	return unique_cycles

#Returns all of the unique pairs of nonrotable bonds in following format (smaller atom #, larger atom #)
def nonrotable_pairs(unique_cycles, bond_matrix):
	pairs = []
	#add all pairs from the cycles
	for cycle in unique_cycles:
		pairs.append([cycle[0], cycle[-1]])
		for i in range(1, len(cycle)):
			pairs.append([cycle[i-1], cycle[i]])

	#add all pairs from the double and tripple bonds
	for i, bonds in enumerate(bond_matrix):
		if 2 in bonds:
			indices_2 = [i for i, x in enumerate(bonds) if x == 2]
			for j in indices_2:
				pairs.append([i + 1, j + 1])
		if 3 in bonds:
			indices_3 = [i for i, x in enumerate(bonds) if x == 2]
			for j in indices_3:
				pairs.append([i + 1, j + 1])

	[x.sort() for x in pairs]
	unique_pairs = []
	for p in pairs:
		if p not in unique_pairs:
			unique_pairs.append(p)
	return unique_pairs

#This function filters out all keys that contain nonrotable pairs in the center of the dihderal
def filtered_key_nonrotable(key, nonrotable_pairs):
	for k in key:
		center = k[1:-1]
		center.sort()
		if center in nonrotable_pairs:
			key.remove(k)
	return key
				
#Generates network based on provided bonds
def generate_network(bond_matrix):
	network = []
	for i in range(0,len(bond_matrix)):
		new_node = Node(i + 1)
		add_neighbors(bond_matrix, new_node)
		network.append(new_node)
	return network

#Adds neightbor to respective node in network
def add_neighbors(bond_matrix, node):
	column = node.atom - 1
	for i in range(0,len(bond_matrix)):
		bond_order = bond_matrix[i][column]
		if bond_order != 0 and bond_order != None:
			node.add_neighbor(i + 1)

		bond_order = bond_matrix[column][i]
		if bond_order != 0 and bond_order != None:
			node.add_neighbor(i + 1)

#Finds all paths of length 4 within the provided network of bonds. 
def len4(bond_network):
	paths = []
	invalid_terminal = []

	for i in range(0,len(bond_network)):
		atom = i + 1
		excludeSet = [atom]

		path = pathfinder(bond_network[i], bond_network, excludeSet, 3, invalid_terminal)
		if path != []:
			for i in path:
				paths.append(i)
		invalid_terminal.append(atom)
	paths.sort()
	return paths

#Recursive pathing from the provided node outward. Terminates after 4 steps
def pathfinder(node, bond_network, excludeSet, n, invalid_terminal):
	if n == 0 and node.atom not in invalid_terminal:
		return [[node.atom]]

	paths = None
	key = []
	for i in node.neighbors:
		if i not in excludeSet:
			excludeSet.append(i)
			paths = pathfinder(bond_network[i-1], bond_network, excludeSet, n-1, invalid_terminal)
			if paths != []:
				for x in paths:
					path = x.append(node.atom)
					key.append(x)
			excludeSet.remove(i)
	return key

#Generates the infomration for the "Full" atom
#Coordiantes, edges, cyclic status
def generate_full_atom(coords, bond_network):
	full_atoms = [None] * len(coords)
	for atom in range(0,len(coords)):
		full_atoms[atom] = [coords[atom], bond_network[atom].neighbors, False]

	return full_atoms

#Returns an array of all the atp,s that are connected to or after the atoms in location 1, 2 in the dihedral_atoms [1,2,3,4]
#Note all atoms are returned in no particular order
##NOTE does not have cycle functionality###
def everythingAfter(dihedral_atoms, full_atoms):
	afterEdges = []

	atom_queue = queue.Queue()
	atom_queue.put(dihedral_atoms[1]) #Starting atom

	visted = []
	visted.append(dihedral_atoms[2])
	visted.append(dihedral_atoms[3]) #Atoms in dihdral that remain stationary

	while not atom_queue.empty():
		atom = atom_queue.get()
		atom_edges = full_atoms[atom - 1][1]

		for a in atom_edges:
			if a not in visted:
				atom_queue.put(a)
				afterEdges.append(a)

		visted.append(atom)
	return afterEdges

#Creates a key that coresponds to the dihedralKey for everythingAfter
def everythingAfterKeyGeneration(dihedralKey, full_atoms):
	everythingAfterKey = [None] * len(dihedralKey)
	for i, key in enumerate(dihedralKey):
		everythingAfterKey[i] = everythingAfter(key, full_atoms)
	return everythingAfterKey

#Creates a dihedral key which links atoms to specified dihedral angles
def dihedral_key_generation(bond_matrix, coords):
	dihedralKey = []
	bond_network = generate_network(bond_matrix)
	full_atoms = generate_full_atom(coords, bond_network)
	dihedralKey = len4(bond_network)

	#Filter out nonrotable dihedral angles
	cycles = cycle_check(full_atoms)
	nonrotable = nonrotable_pairs(cycles, bond_matrix)
	dihedralKey = filtered_key_nonrotable(dihedralKey, nonrotable)

	everythingAfterKey = everythingAfterKeyGeneration(dihedralKey, full_atoms)
	return dihedralKey, everythingAfterKey

#Calcualtes the vector between two atoms
def find_vector(atom1, atom2):
	delta_x = atom1[1] - atom2[1]
	delta_y = atom1[2] - atom2[2]
	delta_z = atom1[3] - atom2[3]
	return [delta_x, delta_y, delta_z]

#Calculates the end vectors for the dihedral angle
def end_vectors(coords, atom_order):
    vector1 = find_vector(coords[atom_order[1] - 1], coords[atom_order[0] - 1])
    vector2 = find_vector(coords[atom_order[2] - 1], coords[atom_order[1] - 1])
    vector3 = find_vector(coords[atom_order[3] - 1], coords[atom_order[2] - 1]) 
    return vector1, vector2, vector3

#Calculates the dihderal angles between the two provided vectors
#https://www.kite.com/python/answers/how-to-get-the-angle-between-two-vectors-in-python
def vector_angle(vector1, vector2):
	unit_vector_1 = vector1 / np.linalg.norm(vector1)
	unit_vector_2 = vector2 / np.linalg.norm(vector2)
	dot_product = np.dot(unit_vector_1, unit_vector_2)
	angle = np.arccos(dot_product)
	return np.degrees(angle)

#Calcualtes all of the dihdral angles for the structures based on the provided key.
def FindDihedral(coords, key):
    individual = [None] * len(key) 
    for i in range(0, len(key)):
        vector1, vector2, vector3 = end_vectors(coords, key[i])
        n1 = np.cross(vector1, vector2)
        n2 = np.cross(vector2, vector3)
        dihedral_angle = vector_angle(n1, n2)
        individual[i] = dihedral_angle
    return individual
