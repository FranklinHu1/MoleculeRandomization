# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 22:19:42 2019

@author: Frank
"""
'''
Algorithm for optimizing the cartesian coordinates of the atoms in the molecule
based on constraints for bond angles and bond lengths. Taking the initial bond angles
and bond lengths, we distort them by adding some random factor. Then, we optimize the 
cartesian coordinates of the molecule using least squares to get as close to the
target set of bond lengths and bond angles as possible, using targetOutputIterate as 
part of the residual calculation function. The resulting "optimized" geometry is 
pickled into a .p file and is used for running gaussian single point later in the 
workflow.
'''

import numpy as np
from scipy.optimize import least_squares
from copy import deepcopy
import random
from targetOutputBase import runGetAnglesAndBondsBase
from targetOutputIterate import getAngleAndBondsIterate
import os
import pickle

# =============================================================================
# Helper Functions and Variables
def cartesianMatrixCorrection(matrix):
    cartesianMatrixNoString = deepcopy(matrix)
    for i in range(len(cartesianMatrixNoString[0])):
        cartesianMatrixNoString[0,i] = float(i)
    return cartesianMatrixNoString

def matrixStrip(matrix):
    matrixResult = np.empty((1, len(matrix)))
    for i in range(len(matrix)):
        matrixResult[0,i] = matrix[i,1]
    return matrixResult

def randomizationOfMatrix(Entry, lowerbound, upperbound):
    MatrixCopy = deepcopy(Entry)
    for row in range(len(MatrixCopy)):
        MatrixCopy[row, 1] += random.uniform(lowerbound, upperbound)
    return MatrixCopy
# =============================================================================
    

def runFunction(startingBondMatrix, startingAngleMatrix, matrices, savePath, time):
    #Could probably move this outside run function to main func, run only once
    workMat = cartesianMatrixCorrection(matrices.originalMatrix)
    flattenedWorkMat = np.array(workMat.flatten(), dtype = object)
    x0 = flattenedWorkMat
    targetAngleMatrix = matrixStrip(randomizationOfMatrix(startingAngleMatrix, -0.1745, 0.1745))
    targetBondMatrix = matrixStrip(randomizationOfMatrix(startingBondMatrix, -0.3, 0.3))
    targetVector = np.concatenate((targetAngleMatrix, targetBondMatrix), axis = None)
    
    def residualCalc(matrix): #Calculates the residuals based off the input matrix
        outPut = getAngleAndBondsIterate(matrix, matrices)
        outAngle, outBond = outPut[0], outPut[1]
        currentAngleMatrix, currentBondMatrix = matrixStrip(outAngle), matrixStrip(outBond)
        currentVector = np.concatenate((currentAngleMatrix, currentBondMatrix), axis = None)
        residuals = np.subtract(targetVector, currentVector)
        return residuals
    
    psol = least_squares(residualCalc, x0, ftol = 1e-3, xtol = 1e-3, gtol = 1e-3)
    
    solution = psol.x.reshape(matrices.numRows, matrices.atomCount)
    solution = np.array(solution, dtype = object)
    for i in range(len(solution[0])):
        for key in matrices.stringToNumDict:
            if matrices.stringToNumDict[key] == float(i):
                solution[0,i] = key
    
    newFilePath = savePath + os.sep + 'resultingGeom' + str(time + 1) + '.p'
    file_out = open(newFilePath, 'wb')
    pickle.dump(solution, file_out)
    file_out.close()

# =============================================================================
# TESTING CODE, ONLY UNCOMMENT WHEN USING
#    beginningAngleMatrix, beginningBondMatrix = matrixStrip(startingAngleMatrix).flatten(), matrixStrip(startingBondMatrix).flatten()
#    beginningVector = np.concatenate((beginningAngleMatrix, beginningBondMatrix), axis = None)
#    
#    mat = psol.x 
#    outPut = getAngleAndBondsIterate(mat, matrices)
#    outAngle, outBond = outPut[0], outPut[1]
#    endingAngleMatrix = matrixStrip(outAngle).flatten()
#    endingBondMatrix = matrixStrip(outBond).flatten()
#    endingVector = np.concatenate((endingAngleMatrix, endingBondMatrix), axis = None)
#    
#    beginningResidual = np.subtract(beginningVector, targetVector)
#    endingResidual = np.subtract(endingVector, targetVector)
#    print('Error at beginning', beginningResidual)
#    print('Error at the end', endingResidual)
# =============================================================================


def mainFuncErrorMin(times, path):
    (startingBondMatrix, startingAngleMatrix, matrices) = runGetAnglesAndBondsBase(path)
    savePath = os.path.join(path, 'resultGeoms')
    for i in range(times):
        runFunction(startingBondMatrix, startingAngleMatrix, matrices, savePath, i)
    return None


        




