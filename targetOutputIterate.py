import math
import numpy as np
import random 
import xlrd
from copy import deepcopy
from matlabProcessing import fullPath, moleculeName

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

molecule, cutoff = fullPath + moleculeName + '.xls', fullPath + 'TabulationsofCutoffBondLengths2.xlsx' 
mat = CreateExcelMat(molecule)
BondLengthCutoffMatrix = CreateExcelMat(cutoff)
originalMatrix = MatrixConversion(mat)
WorkingBondLengthCutoffMatrix = MatrixConversion(BondLengthCutoffMatrix)

stringToNumDict = dict()
for i in range(len(originalMatrix[0])):
    stringToNumDict[originalMatrix[0,i]] = float(i)

#flattenedWorkMat = np.array(([0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,1.8974,4.7545,1.1761,1.1605,1.1761,5.4758,5.4913,5.4757,-0.0055747,-0.013266,1.9321,-0.98217,-0.96077,0.94883,0.95637,-1.951,-0.012094,0.0011149,-0.0030992,1.6547,-1.6981,1.683,-1.6697,0.00011338]),dtype = object)
#template = ['C1C2', 'C1H1', 'C1H2', 'C1H3', 'C2H4', 'C2H5', 'C2H6']

def getAngleAndBondsIterate(matrix, template):
    
#    funcMatrix = deepcopy(matrix)
    funcMatrix = np.array(matrix, dtype = object)
    funcMatrix = funcMatrix.reshape(4 , len(originalMatrix[0])) 
    for i in range(len(funcMatrix[0])):
        for key in stringToNumDict:
            if stringToNumDict[key] == float(i):
                funcMatrix[0,i] = key
    
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
#    
    BondIndexList = template 

    
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

    return [AngleMatrix, BondLengthMatrix]

#print(getAngleAndBondsIterate(flattenedWorkMat, template)) #LINES FOR TESTING















