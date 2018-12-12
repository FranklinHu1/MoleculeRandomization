import math
import numpy as np
from scipy.optimize import least_squares
import random 
import xlrd
from copy import deepcopy
from targetOutputBase import getAngleAndBondsBase, randomizationOfMatrix, originalMatrix
from targetOutputIterate import getAngleAndBondsIterate, stringToNumDict



# =============================================================================
# Helper Functions and Variables
def matrixCorrection(matrix, dictionary):
    workingCopy = deepcopy(matrix)
    matDimensions = workingCopy.shape
    if min(matDimensions) != 2:
        return 'Illegal Matrix Entry'
    for i in range(len(workingCopy)):  
        tempList = []
        for key in dictionary:
            if key in workingCopy[i,0]:
               tempList.append(dictionary[key])
        workingCopy[i,0] = tempList             
    return workingCopy

def cartesianMatrixCorrection(matrix):
    cartesianMatrixNoString = deepcopy(matrix)
    for i in range(len(cartesianMatrixNoString[0])):
        cartesianMatrixNoString[0,i] = float(i)
    return cartesianMatrixNoString

def matrixStrip(matrix):
    matrixCopy = deepcopy(matrix)
    matrixResult = np.empty((1, len(matrixCopy)))
    for i in range(len(matrixCopy)):
        matrixResult[0,i] = matrixCopy[i,1]
    return matrixResult

# =============================================================================
    
#molecule, cutoff = 'Ethane.xlsx','TabulationsofCutoffBondLengths.xlsx'
#startingOutput = getAngleAndBonds(originalMatrix)

def run(originalMatrix = originalMatrix):
    if len(originalMatrix) == 2: #Corrects for molecules that only have either x or only x y coordinates
        multiplier = len(originalMatrix[0])
        originalMatrix = np.append(originalMatrix,[[0]*multiplier],axis = 0)
        originalMatrix = np.append(originalMatrix,[[0]*multiplier],axis = 0)
    elif len(originalMatrix) == 3:
        multiplier = len(originalMatrix[0])
        originalMatrix = np.append(originalMatrix,[[0]*multiplier],axis = 0)
    
    workMat = cartesianMatrixCorrection(originalMatrix)
    #print(workMat)
    flattenedWorkMat = np.array(workMat.flatten(), dtype = object)
    x0 = np.array(flattenedWorkMat, dtype = object)
    #print(flattenedWorkMat)
    startingOutput = getAngleAndBondsBase(flattenedWorkMat) #This line is great for debugging!!
    template = startingOutput[2]
    
    #currentAngleMatrix, currentBondMatrix = matrixStrip(startingOutput[0]), matrixStrip(startingOutput[1])
    #currentCartesianMatrix = startingOutput[2]
    targetAngleMatrix = matrixStrip(randomizationOfMatrix(startingOutput[0], -0.1745, 0.1745)) 
    targetBondMatrix = matrixStrip(randomizationOfMatrix(startingOutput[1],-0.3, 0.3))
    targetVector = np.concatenate((targetAngleMatrix, targetBondMatrix), axis = None) 
    #Ensure that target / current vectors always are 1D of floats, and they have the angle information before the bond information
    
    def residualCalc(matrix): #Insert a counter of some sort to differentiate between the first iteration and later ones
        #In fact, rework the residualCalc functio to get around this dimension error
        mat = deepcopy(matrix)
        outPut = getAngleAndBondsIterate(mat, template)
        rawBondMatrix = outPut[1]
        tempMatrix = []
        for pairs in template:
            for row in range(len(rawBondMatrix)):
                if pairs == rawBondMatrix[row,0]:
                    tempMatrix.append(rawBondMatrix[row])
        tempMatrix = np.array(tempMatrix)
        currentAngleMatrix, currentBondMatrix = matrixStrip(outPut[0]),matrixStrip(tempMatrix)
        currentVector = np.concatenate((currentAngleMatrix, currentBondMatrix), axis = None)
        residuals = np.subtract(targetVector, currentVector)
        return residuals
    
    
    psol = least_squares(residualCalc,x0, ftol = 1e-3, xtol = 1e-3, gtol = 1e-3)
# =============================================================================
#    solution = psol.x.reshape((4 , len(originalMatrix[0]))) 
#    solution = np.array(solution, dtype = object)
#    for i in range(len(solution[0])):
#        for key in stringToNumDict:
#            if stringToNumDict[key] == float(i):
#                solution[0,i] = key
#    
#
#    
#    beginningAngleMatrix, beginningBondMatrix = matrixStrip(startingOutput[0]).flatten(), matrixStrip(startingOutput[1]).flatten()
#    beginningVector = np.concatenate((beginningAngleMatrix, beginningBondMatrix), axis = None)   
# =============================================================================

    
    mat, template = psol.x, template
    outPut = getAngleAndBondsIterate(mat, template)
    rawBondMatrix = outPut[1]
    tempMatrix = []
    for pairs in template:
        for row in range(len(rawBondMatrix)):
            if pairs == rawBondMatrix[row,0]:
                tempMatrix.append(rawBondMatrix[row])
    tempMatrix = np.array(tempMatrix)
    
# =============================================================================
#    endingBondMatrix = matrixStrip(tempMatrix).flatten()
#    endingAngleMatrix = matrixStrip(outPut[0]).flatten()
#    endingVector = np.concatenate((endingAngleMatrix, endingBondMatrix), axis = None)
#    
#    for i in range(len(beginningVector)):
#        beginningVector[i] = abs(beginningVector[i])
#        endingVector[i] = abs(endingVector[i])
#        targetVector[i] = abs(targetVector[i])
#    
#    beginningResidual = np.subtract(beginningVector, targetVector)
#    endingResidual = np.subtract(endingVector, targetVector)
#    print(beginningResidual,'\n', endingResidual)
# =============================================================================


    return [tempMatrix]
#run()






        









#all that remains is to test the result to be sure this is what you want, and then to correct the output

    
    







#sol = least_squares(residualCalc, flattenedWorkMat)

    
#For least squares:
#    1. Calculate the residual vector that takes into account the bond angles and bond lengths
#        Combine the current bond + angles (stripped) and target bond + angles (stripped) into one matrix of matrices
#        Can you subtract the two using np.subtract?
#    2. Pass that vector into the least_squares function
#    3. The initial guess for the least_squares function has to be a vector of floats:
#        To do to this, you can use the functions matrix.flatten() to convert to a single vector
#        You can np.reshape() within the actual function to transform the matrix into a digestable input for your algorithm, it's just to get past the problems with least_squares
    

    
    








    
    










    





    
    
    













## =============================================================================
## MASTER PRINT STATEMENTS DO NOT DELETE
#print('This is the final Bond Matrix with all EXISTING bonds: \n',BondMatrixFinal)
#print('This is the beginning cartesian coordinates of all the atoms (in terms of the bohr radius): \n',workingMatrix)
#print('This is the final Angle matrix with all existing angles: \n',AngleMatrix)    
#print('This is the randomized Angle matrix: \n',AngleMatrixRandomized)
#print('This is the randomized Bond matrix: \n',BondMatrixRandomized)
## =============================================================================
#
#
##'C2H6O.xlsx'
##'TabulationsofCutoffBondLengths.xlsx'