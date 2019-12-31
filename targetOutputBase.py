# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 22:46:54 2019

@author: Frank
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 20:50:33 2019
Beginning initialization for the minimization process
@author: Frank
"""

'''
This module is intended to take the cartesian coordinates of the atoms in a 
molecule, stored in a pickle file called a cartesian matrix, and uses those
cartesian coordinates to calculate all the bond lengths and bond angles in the
molecule. The cartesian coordinates originally come from the minimum energy
geometries of the molecules after running Gaussian on them.
'''
import math
import numpy as np
import xlrd
from copy import deepcopy
import re
import os
import pickle

# =============================================================================
# Helper Functions / variables

def bondLengthDistance(x1,y1,z1,x2,y2,z2):
    return (((x2 - x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))**0.5

BohrRadius = 0.529177 #Given in angstroms

def angleCalculationbyCosineLaw(a,b,c):
    return (((a**2) + (b**2)) - (c**2)) / (2*a*b)

def matrixConversion(excelmat):
    result = np.empty((excelmat.nrows, excelmat.ncols),dtype = object)
    for row in range(excelmat.nrows):
        for col in range(excelmat.ncols):
            result[row,col] = excelmat.cell_value(row,col)
    return result

def createExcelMat(path):
    temporaryMat = xlrd.open_workbook(path)
    result = temporaryMat.sheet_by_index(0)
    return result

def dimensionCorrection(originalMatrix):
     multiplier = len(originalMatrix[0])
     if len(originalMatrix) == 4: return originalMatrix
     elif len(originalMatrix) == 2:
         originalMatrix = np.append(originalMatrix, [[0] * multiplier], axis = 0)
         originalMatrix = np.append(originalMatrix, [[0] * multiplier], axis = 0)
     elif len(originalMatrix) == 3:
         originalMatrix = np.append(originalMatrix, [[0] * multiplier], axis = 0)
     return originalMatrix
             
# =============================================================================

class moleculeAndVariables: #Object that initalizes the constants for this molecule
    def __init__(self, path):
        self.originalMatrixPath = path + os.sep + 'moleculeCartMatrix.p'
        objFile = open(self.originalMatrixPath, 'rb')
        matrix = pickle.load(objFile)
        self.originalMatrix = np.array(matrix, dtype = object)
        self.cutoff = os.getcwd() + os.sep + 'TabulationsofCutoffBondLengths2.xlsx'
        self.originalMatrix = dimensionCorrection(self.originalMatrix)
        bondLengthCutoffMatrix = createExcelMat(self.cutoff)
        self.workingBondLengthCutoffMatrix = matrixConversion(bondLengthCutoffMatrix)
        self.atomCount = len(self.originalMatrix[0])
        self.numRows = 4
    
        
def initializeMatrices(path): #Called by main func to initalize the matrices
    matrices = moleculeAndVariables(path)
    return matrices

def createStringToNumDict(matrix): #Takes the original matrix and creates the string to num dictionary
    stringToNumDict = dict()
    for i in range(len(matrix[0])):
        stringToNumDict[matrix[0,i]] = float(i)
    return stringToNumDict

def createBondLengthCutoffDict(cutoffMat): #Creates a dictionary of the bondlengthcutoffs
    bondLengthCutoffs = dict() 
    for row in range(1, len(cutoffMat)):
        key , value = cutoffMat[row,0] , cutoffMat[row,1:]
        bondLengthCutoffs[key] = value
    return bondLengthCutoffs


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


def addInBondLengths(bondLengthMatrix, dimensions, funcMatrix): #Destructively adds things into the bond length function
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


def removeBadRows(bondLengthMatrix, bondLengthCutoffs): #THIS IS THE BLOCK THAT SORTS THROUGH BOND LENGTHS
    badRows, numberBadRows = [], 0
    for row in range(len(bondLengthMatrix)): #For all bonds that exceed the bond length, the value is set to None
        Checkmatrix = []
        for character in bondLengthMatrix[row,0]:
            if character.isalpha():
                Checkmatrix.append(character)
        Checkstring = Checkmatrix[0] + Checkmatrix[1]
        Threshold = bondLengthCutoffs[Checkstring]
        if (math.isclose(bondLengthMatrix[row, 1], Threshold[0], abs_tol = 1e-1) == False) and\
        (math.isclose(bondLengthMatrix[row, 1], Threshold[1], abs_tol = 1e-1) == False) and\
        (math.isclose(bondLengthMatrix[row, 1], Threshold[2], abs_tol = 1e-1) == False):
            badRows.append(row)
            numberBadRows += 1
    numberGoodRows = bondLengthMatrix.shape[0] - numberBadRows #Takes the number of total rows and subtracts out the number of bad rows from the matrix
    bondLengthMatrixFinal = np.empty((numberGoodRows, 2), dtype = object)
    rowCount = 0
    for row in range(len(bondLengthMatrix)): #Brings all legal bonds into a final matrix
        if row not in badRows:
            bondLengthMatrixFinal[rowCount, 0], bondLengthMatrixFinal[rowCount, 1] = bondLengthMatrix[row, 0], bondLengthMatrix[row, 1]
            rowCount += 1
    return bondLengthMatrixFinal
        

def findMultiBondedAtoms(bondLengthMatrix):
    bondIndexList = bondLengthMatrix[:, 0].tolist()
    uniqueAtomSet, allAtomList = set(), []
    for i in range(len(bondIndexList)):
        totalcharacterlist = list(bondIndexList[i])
        symbolIndex = []
        for index in range(len(totalcharacterlist)):
            if totalcharacterlist[index].isalpha():
                symbolIndex.append(index)
        Atom1, Atom2 = '',''
        for j in range(symbolIndex[1]):
            Atom1 += totalcharacterlist[j]
        for h in range(symbolIndex[1], len(totalcharacterlist)):
            Atom2 += totalcharacterlist[h]
        uniqueAtomSet.add(Atom1)
        uniqueAtomSet.add(Atom2)
        allAtomList.append(Atom1)
        allAtomList.append(Atom2)
    multiBondedAtoms = []
    for s in uniqueAtomSet:
        count = allAtomList.count(s)
        if count > 1:
            multiBondedAtoms.append(s)
    return multiBondedAtoms, bondIndexList
    

def figureAtomTriplets(multiBondedAtomsList, bondIndexList):
    tripleList = [] #TripleList is a list, so you have to use [row][col] rather than the nice [row,col] you get with a np array
    for atom in multiBondedAtomsList:
        Indlist = []
        for index in range(len(bondIndexList)):
            if atom in bondIndexList[index]:
                Indlist.append(index)
        for i in range(len(Indlist)):
            for j in range(i+1, len(Indlist)):
                Atmstr1 = bondIndexList[Indlist[i]].split(atom)
                Atmstr2 = bondIndexList[Indlist[j]].split(atom)
                Temp1, Temp2 = None, None
                for IND in range(len(Atmstr1)):
                    if Atmstr1[IND] != '':
                        Temp1 = Atmstr1[IND]
                    if Atmstr2[IND] != '':
                        Temp2 = Atmstr2[IND]
                triplet = [Temp1, atom, Temp2]
                tripleList.append(triplet)
    return tripleList


def testInTripleList(bondMatrixEntry, tripleListEntry):
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
    

def calculateAngles(bondMatrixRaw, tripleList, funcMatrix):
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

def getAngleAndBondsBase(matrices, stringToNumDict):
    funcMatrix = np.array(matrices.originalMatrix, dtype = object)
    for i in range(len(funcMatrix[0])): 
        for key in stringToNumDict:
            if stringToNumDict[key] == float(i):
                funcMatrix[0,i] = key
    bondLengthCutoffs = createBondLengthCutoffDict(matrices.workingBondLengthCutoffMatrix)
    bondMatrixReturns = createEmptyBondMatrix(funcMatrix)
    bondMatrix, workingMatrixDimensions = bondMatrixReturns[0], bondMatrixReturns[1]
    addInBondLengths(bondMatrix, workingMatrixDimensions, funcMatrix)
    bondMatrixRaw = deepcopy(bondMatrix)
    bondMatrixFinal = removeBadRows(bondMatrix, bondLengthCutoffs)
    multiBondedAtomsResult = findMultiBondedAtoms(bondMatrixFinal)
    multiBondedAtoms, bondIndexList = multiBondedAtomsResult[0], multiBondedAtomsResult[1]
    tripleList = figureAtomTriplets(multiBondedAtoms, bondIndexList)
    angleMatrix = calculateAngles(bondMatrixRaw, tripleList, funcMatrix)
    bondLengthTemplate = bondMatrixFinal[:, 0].tolist()
    return bondMatrixFinal, angleMatrix, bondLengthTemplate, tripleList


def runGetAnglesAndBondsBase(path):
    matrices = initializeMatrices(path)
    stringToNumDict = createStringToNumDict(matrices.originalMatrix)
    matrices.stringToNumDict = stringToNumDict
    output = getAngleAndBondsBase(matrices, stringToNumDict)
    bondMatrixFinal, angleMatrix = output[0], output[1]
    bondLengthTemplate, tripleList = output[2], output[3]
    matrices.bondLengthTemplate = bondLengthTemplate
    matrices.angleTemplate = tripleList
    return bondMatrixFinal, angleMatrix, matrices

    
    
                
