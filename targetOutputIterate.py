# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 10:35:47 2019
Function called within the residual calc function for error minimization
@author: Frank
"""
'''
This is a modified version of targetOutputBase.py, and it's used in the iterative
minimization of errors. It still takes in the parameters (the cartesians of the 
atoms in the molecule) and returns the bond lengths and bond angles from 
calculations performed on those cartesians.
'''

import math
import numpy as np
import re

# =============================================================================
# Helper Functions / variables

def bondLengthDistance(x1,y1,z1,x2,y2,z2):
    return (((x2 - x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))**0.5

BohrRadius = 0.529177 #Given in angstroms

def angleCalculationbyCosineLaw(a,b,c):
    return (((a**2) + (b**2)) - (c**2)) / (2*a*b)

# =============================================================================


def addInBondValuesBasedOnTemplate(bondLengthMatrix, template, funcMatrix): #Takes in the template and creates the final bond matrix
    place = 0
    for bonds in template:
        tokenList = re.split('(\d+)', bonds)
        atom1 = tokenList[0] + tokenList[1]
        atom2 = tokenList[2] + tokenList[3]
        usefulCols = []
        for i in range(len(funcMatrix[0])):
            if atom1 == funcMatrix[0, i] or atom2 == funcMatrix[0, i]:
                usefulCols.append(funcMatrix[1:, i].tolist())
        X1, X2 = usefulCols[0][0], usefulCols[1][0]
        Y1, Y2 = usefulCols[0][1], usefulCols[1][1]
        Z1, Z2 = usefulCols[0][2], usefulCols[1][2]
        bondLengthNotAngstrom = bondLengthDistance(X1, Y1, Z1, X2, Y2, Z2)
        bondLength = bondLengthNotAngstrom 
        bondLengthMatrix[place, 0] = bonds
        bondLengthMatrix[place, 1] = bondLength
        place += 1

def addInBondLengthsRaw(bondLengthMatrix, dimensions, funcMatrix): #Destructively adds things into the bond length function
    place = 0
    for i in range(dimensions[1]):
        if 'H' in funcMatrix[0,i]:
            pass
        else:
            for j in range(i + 1, dimensions[1]):
                index = funcMatrix[0,i] + funcMatrix[0,j]
                bondLengthMatrix[place,0] = index
                X1,X2 = funcMatrix[1,i], funcMatrix[1,j]
                Y1,Y2 = funcMatrix[2,i], funcMatrix[2,j]
                Z1,Z2 = funcMatrix[3,i], funcMatrix[3,j] 
                BondLength = bondLengthDistance(X1,Y1,Z1,X2,Y2,Z2)
                bondLengthMatrix[place,1] = BondLength 
                place += 1
    pass

def createEmptyBondMatrix(matrix): #Creates the empty bond matrix of specified number of rows
    workingMatrix = matrix
    workingMatrixDimensions = workingMatrix.shape
    headingList = workingMatrix[0,:].tolist()
    atomString = ''
    for i in range(len(headingList)):
        atomString += headingList[i]
    numberHydrogens = atomString.count('H') #Gives the # of Hydorgens in the molecule, which necessarily cannot form bonds with each other
    numberNonHydrogens = workingMatrixDimensions[1] - numberHydrogens #workingMatrixDimensions[1] is equivalent to the total number of atoms
    bondLengthMatrixRows = 0
    #You start with the maximum number of bonds (the first atom can bond to all atoms after it) and work down by one each time. For example,
    #   in ethane, you have C2H6. The first carbon can form 7, and the second carbon can form 6 (we already counted C1C2!), and we continue
    #   adding by one less until we finish with the number of non-hydrogen atoms
    for i in range(workingMatrixDimensions[1] - 1, (workingMatrixDimensions[1]-numberNonHydrogens)-1, -1): 
        bondLengthMatrixRows += i
    bondLengthMatrix = np.empty((bondLengthMatrixRows, 2), dtype = object)
    return bondLengthMatrix, workingMatrixDimensions


def testInTripleList(bondMatrixEntry, tripleListEntry): #Sees if a bond can be found in a triplet
    splitCharacters = re.split('(\d+)', bondMatrixEntry)
    atom1 = splitCharacters[0] + splitCharacters[1]
    atom2 = splitCharacters[2] + splitCharacters[3]
    if atom1 in tripleListEntry and atom2 in tripleListEntry: return True
    else: return False

def addingSides(tripleList, index, funcMatrix):
    needDistanceAtoms = []
    needDistanceAtoms.append(tripleList[index][0])
    needDistanceAtoms.append(tripleList[index][2])
    usefulCol = []
    for i in range(len(needDistanceAtoms)):
        search = needDistanceAtoms[i]
        for col in range(len(funcMatrix[0])):
            if funcMatrix[0, col] == search: #We found the column of one of the side atoms
                usefulCol.append(col)
    X1,X2 = funcMatrix[1, usefulCol[0]], funcMatrix[1, usefulCol[1]]
    Y1,Y2 = funcMatrix[2, usefulCol[0]], funcMatrix[2, usefulCol[1]]
    Z1,Z2 = funcMatrix[3, usefulCol[0]], funcMatrix[3, usefulCol[1]]
    bondLength = bondLengthDistance(X1, Y1, Z1, X2, Y2, Z2)
    key, value = (needDistanceAtoms[0] + needDistanceAtoms[1]), bondLength 
    return key, value


def calculateAngles(bondMatrixRaw, tripleList, funcMatrix): #Function for calculating the angles
    angleMatrix = np.empty((len(tripleList),2),dtype = object)
    for i in range(len(tripleList)):
        sides = dict()
        for j in range(len(bondMatrixRaw)): #We use the raw bond matrix because we want all the sides
            bondMatrixEntry = bondMatrixRaw[j, 0]
            tripleListEntry = tripleList[i]
            if testInTripleList(bondMatrixEntry, tripleListEntry):
                key, value = bondMatrixRaw[j, 0], bondMatrixRaw[j, 1]
                sides[key] = value
        if len(sides) != 3:
            addingSidesResult = addingSides(tripleList, i, funcMatrix)
            key, value = addingSidesResult[0], addingSidesResult[1]
            sides[key] = value
        CentralAtom = tripleList[i][1]
        NonOpposingSides, c = [], None #the variable c represents the hypotenuse in a law of cosines calculation
        for key in sides:
            if CentralAtom not in key:
                c = sides[key]
            else:
                NonOpposingSides.append(sides[key])
        a = NonOpposingSides[0]
        b = NonOpposingSides[1]
        UnconvertedAngle = angleCalculationbyCosineLaw(a,b,c)
        if math.isclose(UnconvertedAngle, -1, abs_tol = 1e-2): #Handles the case when the angle is approximately pi radians
            Angle = math.pi
        elif math.isclose(UnconvertedAngle, 1, abs_tol = 1e-2): #Handles the case (if it happens) where the angle is 0 radians
            Angle = 0
        else:
            Angle = math.acos(UnconvertedAngle)
        angleMatrix[i,0], angleMatrix[i,1]  = tripleList[i], Angle
    return angleMatrix


# =============================================================================
# TESTING CODE, comment out when not using

#template = ['C1C2', 'C1H1', 'C1H2', 'C1H3', 'C2H4', 'C2H5', 'C2H6']
#bondLengthMatrix = np.empty((len(template), 2), dtype = object)
#funcMatrix = flattenedWorkMat.reshape(4, 8)
#funcMatrix[0, :] = ['C1', 'C2', 'H1', 'H2', 'H3', 'H4', 'H5', 'H6']
#dimensions = funcMatrix.shape
#
#emptyBondMatrix = createEmptyBondMatrix(funcMatrix)[0]
#addInBondLengthsRaw(emptyBondMatrix, dimensions, funcMatrix)
#tripleList = [['C1', 'C2', 'H4'], ['C1', 'C2', 'H5'], ['C1', 'C2', 'H6'], ['H4', 'C2', 'H5'], ['H4', 'C2', 'H6'], ['H5', 'C2', 'H6'], ['C2', 'C1', 'H1'], ['C2', 'C1', 'H2'], ['C2', 'C1', 'H3'], ['H1', 'C1', 'H2'], ['H1', 'C1', 'H3'], ['H2', 'C1', 'H3']]
#print(addInBondValuesBasedOnTemplate(bondLengthMatrix, template, funcMatrix))

# =============================================================================


def getAngleAndBondsIterate(matrix, matrices): #Take the matrices object from targetOutputBaseNew.py
    funcMatrix = np.array(matrix, dtype = object)
    funcMatrix = funcMatrix.reshape(matrices.numRows, matrices.atomCount)
    for i in range(len(funcMatrix[0])): 
        for key in matrices.stringToNumDict:
            if matrices.stringToNumDict[key] == float(i):
                funcMatrix[0,i] = key
    bondLengthMatrixFinal = np.empty((len(matrices.bondLengthTemplate), 2), dtype = object)
    rawBondMatrxiResult = createEmptyBondMatrix(funcMatrix)
    bondMatrixRaw, dimensions = rawBondMatrxiResult[0], rawBondMatrxiResult[1]
    addInBondValuesBasedOnTemplate(bondLengthMatrixFinal, matrices.bondLengthTemplate, funcMatrix)
    addInBondLengthsRaw(bondMatrixRaw, dimensions, funcMatrix)
    tripleList = matrices.angleTemplate
    angleMatrix = calculateAngles(bondMatrixRaw, tripleList, funcMatrix)
    return angleMatrix, bondLengthMatrixFinal
    
    

    
    
    
    
    
        
    
    
    
    
    