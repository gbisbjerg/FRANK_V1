import numpy as np
import copy


#Calculates the euclidean distance between two points
def atomDistance(atom_a, atom_b):
    return np.sqrt( (atom_a[0] - atom_b[0])**2 + (atom_a[1] - atom_b[1])**2 + (atom_a[2] - atom_b[2])**2 )

#Prints constrain matrix 
def print_Contraint_Matrix(constraintMatrix):
    print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in constraintMatrix]))

#Generates a constrain matrix based on atom distances, pairs will be evaluated as follows
    #  0 No constrain violations
   # -1 Atoms are too close
    #  1 Atoms are too far away
def RDGA_Contraint_Matrix(coordinates, bonds):
    noBondMinDistance = 1.5
    BondMaxDistance = 2.2
    BondMinDistance = 0.95

    numAtoms = len(coordinates) #Length of 
    #Initialize matrix to track contraints
    contraintViolations = [None] * numAtoms

    for i in range(numAtoms): #Rows
        atomRow = [None] * numAtoms
        for j in range(numAtoms): #Columns
            if i <= j: #Sets all pairs on or above the diagonal to zero
                atomRow[j] = 0
            else: #if the atoms are different
                distance = atomDistance(coordinates[i], coordinates[j])
                if bonds[i][j] != 0: #If there is a bond between the two atoms
                    if distance > BondMaxDistance:
                        atomRow[j] = 1
                    elif distance < BondMinDistance:
                        atomRow[j] = -1
                    else:
                        atomRow[j] = 0
                else: #if there is no bond located between the two atoms
                    if distance < noBondMinDistance: 
                        atomRow[j] = -1
                    else:
                        atomRow[j] = 0
        contraintViolations[i] = atomRow
    return contraintViolations

#Function to translate coordiantes of an atom pair
    #atoms are moved closer if the step is postive
    #futher away if the step is negative
def translate(bond, coordinates, step):
    #Note the step amount is split between both atoms in the bond
    half_step = step/2

    #vector is noted by bond[0] -> bond[1]
    vector = np.subtract(coordinates[bond[0]], coordinates[bond[1]])
    unit_vector = np.divide(vector, np.linalg.norm(vector))
    step_vector = np.multiply(unit_vector,  half_step)

    #Translate atom0 (bond[0]) in direction of tail and atom1 in direction of head of the vector
    coordinates[bond[0]] = np.add(coordinates[bond[0]], np.multiply(step_vector, -1)).tolist()
    coordinates[bond[1]] = np.add(coordinates[bond[1]], step_vector).tolist()
    return coordinates

#Convert from xyzCoords back to full coordinates (including atom type)
def xyzCoords_to_gamess(xyzCoords, gamessCoords):
    new_gamessCoords = copy.deepcopy(gamessCoords)
    for i, row in enumerate(new_gamessCoords):
        new_gamessCoords[i] = row[:2] + xyzCoords[i]
    return new_gamessCoords

#Uses the constraintMatrix to translate violations
def RGDA_one_iteration(constraintMatrix, coordinates, step_size):
    for i, row in enumerate(constraintMatrix):
        for j, value in enumerate(row):
            if value == 0:
                pass
            elif value == -1: #Values in bond are too close
                coordinates = translate([i,j], coordinates, -step_size)
            else:
                coordinates = translate([i,j], coordinates, step_size)
    return coordinates

#Evaluates if there are any constraints in the constraintMatrix that are violated
def ConstrainCheck(constraintMatrix):
    for row in constraintMatrix:
        np_row = np.array(row)
        is_all_zero = np.all((np_row == 0))
        if not is_all_zero:
            return False
    return True

#Runs the RGDA on a single structure
def RGDA_one_structure(xyzCoords, bondMatrix, step_size, iterations):
    #Generate the constraintMatrix

    while iterations > 0: #Driving loop
        constraintMatrix = RDGA_Contraint_Matrix(xyzCoords, bondMatrix)

        if iterations == 1:
            pass
            # print("Attempt ran out of iterations")
            # print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in constraintMatrix]))

        if ConstrainCheck(constraintMatrix):
            # print("Ending on iteration:", iterations)
            # print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in constraintMatrix]))
            break

        xyzCoords = RGDA_one_iteration(constraintMatrix, xyzCoords, step_size)
        iterations -= 1
    return xyzCoords

