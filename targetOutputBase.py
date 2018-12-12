#This version of the getAngleAndBonds program work with a 1D np array of floats

import math
import numpy as np
import random 
import xlrd
from copy import deepcopy
from matlabProcessing import fullPath


# =============================================================================
# Helper Functions / variables

def BondLengthDistance(x1,y1,z1,x2,y2,z2):
    return (((x2 - x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))**0.5

BohrRadius = 0.529177 #Given in angstroms

def AngleCalculationbyCosineLaw(a,b,c):
    return (((a**2) + (b**2)) - (c**2)) / (2*a*b)

def MatrixConversion(excelmat):
    result = np.empty((excelmat.nrows, excelmat.ncols),dtype = object)
    for row in range(excelmat.nrows):
        for col in range(excelmat.ncols):
            result[row,col] = excelmat.cell_value(row,col)
    return result

def CreateExcelMat(path):
    temporaryMat = xlrd.open_workbook(path)
    result = temporaryMat.sheet_by_index(0)
    return result

def randomizationOfMatrix(Entry, lowerbound, upperbound):
    MatrixCopy = deepcopy(Entry)
    for row in range(len(MatrixCopy)):
        MatrixCopy[row, 1] += random.uniform(lowerbound, upperbound)
    return MatrixCopy
# =============================================================================
molecule, cutoff = fullPath + 'C6N1H11.xls', fullPath + 'TabulationsofCutoffBondLengths2.xlsx' #Input a new molecule each time you want to use a different one
mat = CreateExcelMat(molecule)
originalMatrix = MatrixConversion(mat)
BondLengthCutoffMatrix = CreateExcelMat(cutoff)
WorkingBondLengthCutoffMatrix = MatrixConversion(BondLengthCutoffMatrix)


#flattenedWorkMat = np.array(([0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,1.8974,4.7545,1.1761,1.1605,1.1761,5.4758,5.4913,5.4757,-0.0055747,-0.013266,1.9321,-0.98217,-0.96077,0.94883,0.95637,-1.951,-0.012094,0.0011149,-0.0030992,1.6547,-1.6981,1.683,-1.6697,0.00011338]),dtype = object)
#lINE 46 IS FOR TESTING
#Molecule of Ethane, rehsape 4,8


stringToNumDict = dict()
for i in range(len(originalMatrix[0])):
    stringToNumDict[originalMatrix[0,i]] = float(i)

#stringToNumDict[key] == funcMatrix[0,i]
def getAngleAndBondsBase(matrix, cutoffMat = WorkingBondLengthCutoffMatrix):
    
    funcMatrix = deepcopy(matrix)
    funcMatrix = np.array(matrix, dtype = object)
    funcMatrix = funcMatrix.reshape(4 , 18) #Take care to check this line every time you are testing, so that it's corresponding to the correct molecule
    for i in range(len(funcMatrix[0])):
        for key in stringToNumDict:
            if stringToNumDict[key] == float(i):
                funcMatrix[0,i] = key
#            if type(funcMatrix[0,i]) == str:
#                pass
#            else:
#                if math.isclose(stringToNumDict[key],funcMatrix[0,i], rel_tol = 1e-01, abs_tol = 1e-01): #Feel free to adjust the tolerance as necessary
#                    funcMatrix[0,i] = key #There's a problem with this line, I don't know why?
                #Apparently, this funcMatrix does not allow for elements to be reassigned as sr
#    print(funcMatrix)
            
    BondLengthCutoffs = dict() #Creates a dictionary of the bondlengthcutoffs
    for row in range(1, len(cutoffMat)):
        Key , Value = cutoffMat[row,0] , cutoffMat[row,1:]
        BondLengthCutoffs[Key] = Value
    
#    print(BondLengthCutoffs)
    
    workingMatrix = funcMatrix #Important! if you change anything else, look at this line first
    workingMatrixDimensions = workingMatrix.shape
    HeadingList = workingMatrix[0,:]
    AtomString = ''
    for i in range (len(HeadingList)):
        AtomString += HeadingList[i]
        
    NumberHydrogens = AtomString.count('H') #Gives the # of Hydorgens in the molecule, which necessarily cannot form bonds with each other
    NumberNonHydrogens = workingMatrixDimensions[1] - NumberHydrogens #workingMatrixDimensions[1] is equivalent to the total number of atoms
    BondLengthMatrixRows = 0
    
    for i in range (workingMatrixDimensions[1]-1, (workingMatrixDimensions[1]-NumberNonHydrogens)-1, -1):
        BondLengthMatrixRows += i
    BondLengthMatrix = np.empty((BondLengthMatrixRows, 2), dtype = object)
    
    place = 0
    for i in range(workingMatrixDimensions[1]):
        if 'H' in workingMatrix[0,i]:
            pass
        else:
            for j in range(i+1, workingMatrixDimensions[1]):
                index = workingMatrix[0,i] + workingMatrix[0,j]
                BondLengthMatrix[place,0] = index
                X1,X2 = workingMatrix[1,i], workingMatrix[1,j]
                Y1,Y2 = workingMatrix[2,i], workingMatrix[2,j]
                Z1,Z2 = workingMatrix[3,i], workingMatrix[3,j] 
                BondLength = BondLengthDistance(X1,Y1,Z1,X2,Y2,Z2)
                BondLengthMatrix[place,1] = BondLength*BohrRadius
                place += 1
    
    
    NumberBadRows = 0 #THIS IS THE BLOCK THAT SORTS THROUGH BOND LENGTHS
    BadRows = []
    for row in range(len(BondLengthMatrix)): #For all bonds that exceed the bond length, the value is set to None
        Checkmatrix = []
        for character in BondLengthMatrix[row,0]:
            if character.isalpha():
                Checkmatrix.append(character)
        Checkstring = Checkmatrix[0] + Checkmatrix[1]
        Threshold = BondLengthCutoffs[Checkstring]
        if (math.isclose(BondLengthMatrix[row, 1], Threshold[0], abs_tol = 1e-1) == False) and (math.isclose(BondLengthMatrix[row, 1], Threshold[1], abs_tol = 1e-1) == False) and (math.isclose(BondLengthMatrix[row, 1], Threshold[2], abs_tol = 1e-1) == False):
            BadRows.append(row)
    #        BondLengthMatrix[row,1] = None
            NumberBadRows += 1
    
    NumberGoodRows = BondLengthMatrix.shape[0] - NumberBadRows #Takes the number of total rows and subtracts out the number of bad rows from the matrix
    BondMatrixFinal = np.empty((NumberGoodRows, 2), dtype = object)
    
    rowCount = 0
    for row in range(len(BondLengthMatrix)): #Brings all legal bonds into a final matrix
        if row not in BadRows:
            BondMatrixFinal[rowCount,0], BondMatrixFinal[rowCount,1] = BondLengthMatrix[row,0],BondLengthMatrix[row,1]
            rowCount += 1
    
    BondIndexList = [] #THIS IS THE LINE OF ISSUE, FIGURE OUT HOW TO MAKE THIS WORK WITH TEMPLATE!!!!
    #The answer will most likely be in changing the error min subroutine file
    for row in range(len(BondMatrixFinal)):
        BondIndexList.append(BondMatrixFinal[row,0])
    
    UniqueAtomSet = set()
    AllAtomlist = []
    for i in range(len(BondIndexList)):
        totalcharacterlist = list(BondIndexList[i])
        symbolIndex = []
        for index in range(len(totalcharacterlist)):
            if totalcharacterlist[index].isalpha():
                symbolIndex.append(index)
        Atom1, Atom2 = '',''
        for j in range(symbolIndex[1]):
            Atom1 += totalcharacterlist[j]
        for h in range(symbolIndex[1], len(totalcharacterlist)):
            Atom2 += totalcharacterlist[h]
        UniqueAtomSet.add(Atom1)
        UniqueAtomSet.add(Atom2)
        AllAtomlist.append(Atom1)
        AllAtomlist.append(Atom2)
        
    MultiBondedAtoms = [] #These are the atoms that form multiple bonds and can be used as the center atom for calculating bond angles
    for s in UniqueAtomSet:
        count = AllAtomlist.count(s)
        if count > 1:
            MultiBondedAtoms.append(s)
    
    #Creating the triples of bonded atoms!
    #print(BondIndexList)
    TripleList = [] #TripleList is a list, so you have to use [row][col] rather than the nice [row,col] you get with a np array
    for atom in MultiBondedAtoms:
        Indlist = []
        for index in range(len(BondIndexList)):
            if atom in BondIndexList[index]:
                Indlist.append(index)
        for i in range(len(Indlist)):
            for j in range(i+1, len(Indlist)):
                Atmstr1 = BondIndexList[Indlist[i]].split(atom)
                Atmstr2 = BondIndexList[Indlist[j]].split(atom)
                Temp1, Temp2 = None, None
                for IND in range(len(Atmstr1)):
                    if Atmstr1[IND] != '':
                        Temp1 = Atmstr1[IND]
                    if Atmstr2[IND] != '':
                        Temp2 = Atmstr2[IND]
                triplet = [Temp1, atom, Temp2]
                TripleList.append(triplet)
    #Now that you have the matrix with the bond triples, find a way to calculate the actual angle!
    AngleMatrix = np.empty((len(TripleList),2),dtype = object)
    for i in range(len(TripleList)):
        sides = dict()
        for j in range(len(BondLengthMatrix)):
            if (BondLengthMatrix[j,0] == TripleList[i][1] + TripleList[i][0]) or (BondLengthMatrix[j,0] == TripleList[i][0] + TripleList[i][1]):
                key, value = BondLengthMatrix[j,0],BondLengthMatrix[j,1]
                sides[key] = value
            if (BondLengthMatrix[j,0] == TripleList[i][1] + TripleList[i][2]) or (BondLengthMatrix[j,0] == TripleList[i][2] + TripleList[i][1]):
                key, value = BondLengthMatrix[j,0],BondLengthMatrix[j,1]
                sides[key] = value
            if (BondLengthMatrix[j,0] == TripleList[i][0] + TripleList[i][2]) or (BondLengthMatrix[j,0] == TripleList[i][2] + TripleList[i][0]):
                key, value = BondLengthMatrix[j,0],BondLengthMatrix[j,1]
                sides[key] = value
        if len(sides) != 3:
            NeedDistanceAtoms = []
            NeedDistanceAtoms.append(TripleList[i][0])
            NeedDistanceAtoms.append(TripleList[i][2])
            UsefulCol = []
            for IND in range(len(NeedDistanceAtoms)):
                Search = NeedDistanceAtoms[IND]
                for col in range(len(workingMatrix[0])):
                    if workingMatrix[0,col] == Search:
                        UsefulCol.append(col)
            X1,X2 = workingMatrix[1,UsefulCol[0]], workingMatrix[1,UsefulCol[1]]
            Y1,Y2 = workingMatrix[2,UsefulCol[0]], workingMatrix[2,UsefulCol[1]]
            Z1,Z2 = workingMatrix[3,UsefulCol[0]], workingMatrix[3,UsefulCol[1]]
            BondLength = BondLengthDistance(X1,Y1,Z1,X2,Y2,Z2)
            key, value = (NeedDistanceAtoms[0] + NeedDistanceAtoms[1]), BondLength*BohrRadius
            sides[key] = value
        CentralAtom = TripleList[i][1]
        NonOpposingSides,c = [],None
        for key in sides:
            if CentralAtom not in key:
                c = sides[key]
            else:
                NonOpposingSides.append(sides[key])
        a = NonOpposingSides[0]
        b = NonOpposingSides[1]
        UnconvertedAngle = AngleCalculationbyCosineLaw(a,b,c)
        if math.isclose(UnconvertedAngle, -1, abs_tol = 1e-2): #Handles the case when the angle is approximately pi radians
            Angle = math.pi
        elif math.isclose(UnconvertedAngle, 1, abs_tol = 1e-2): #Handles the case (if it happens) where the angle is 0 radians
            Angle = 0
        else:
            Angle = math.acos(UnconvertedAngle)
        AngleMatrix[i,0],AngleMatrix[i,1]  = TripleList[i], Angle
    #Takes the angle matrix and randomizes it within the range given
#    AngleMatrixRandomized = RandomizationofMatrix(AngleMatrix, -0.0872665, 0.0872665)
#    BondMatrixRandomized = RandomizationofMatrix(BondMatrixFinal, -0.3, 0.3)
    
    bondLengthTemplate = []
    for i in range(len(BondMatrixFinal)):
        if BondMatrixFinal[i,0] not in bondLengthTemplate:
            bondLengthTemplate.append(BondMatrixFinal[i,0])
            
    return [AngleMatrix, BondMatrixFinal, bondLengthTemplate]


#print(getAngleAndBondsBase(flattenedWorkMat)) #LINE FOR TESTING


    
    











