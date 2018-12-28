'''
code helps visualize the differences between the original molecule (in blue) against the new molecule (in red)
'''
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from errorMinSubRoutine import run
from targetOutputBase import originalMatrix
import numpy as np
import re


class constant:
    def __init__(self, template, solution, original): #Method of storing the template, solution, and original
        self.template = template #The template of all bonds
        self.solution = solution #The solution pulled from the errorMinSubRoutine
        self.original = original #The original cartesian matrix
        
def originalMatrixCheck(originalMatrix):
    if len(originalMatrix) == 2: #Corrects for molecules that only have either x or only x y coordinates
        multiplier = len(originalMatrix[0])
        originalMatrix = np.append(originalMatrix,[[0]*multiplier],axis = 0)
        originalMatrix = np.append(originalMatrix,[[0]*multiplier],axis = 0)
    elif len(originalMatrix) == 3:
        multiplier = len(originalMatrix[0])
        originalMatrix = np.append(originalMatrix,[[0]*multiplier],axis = 0)
    return originalMatrix
    

def getTemplateAndSolution():
    matrix, solution = run()[0], run()[1]
#    print(matrix)
    template = matrix[:,0]
#    print(template)
    result = originalMatrixCheck(originalMatrix)
    constants = constant(template, solution, result)
    return constants


def figureBondedAtoms(template, origMatrix):
    tuples = []
    for i in range(len(template)):
        string = template[i]
        lst = re.split('(\d+)',string)
        firstAtom, secondAtom = lst[0] + lst[1], lst[2] + lst[3]
        tempList = []
        for j in range(len(origMatrix[0])):
            if origMatrix[0,j] == firstAtom or origMatrix[0,j] == secondAtom:
                tempList.append(j)
        tuples.append(tuple(tempList))
#    print(tuples)
    return tuples

def gatherCartesianData(tupleList, origMatrix):
    cartesianPlottingMatrix = []
    for i in range(len(tupleList)):
        firstAtomIndex, secondAtomIndex = tupleList[i][0], tupleList[i][1]
        firstAtomCoordinates = np.array(origMatrix[1:,firstAtomIndex], dtype = float)
        secondAtomCoordinates = np.array(origMatrix[1:,secondAtomIndex], dtype = float)
        dataTuple = (firstAtomCoordinates, secondAtomCoordinates)
        cartesianPlottingMatrix.append(dataTuple)
#    print(cartesianPlottingMatrix)
    return cartesianPlottingMatrix

def plotting(matrixOriginal, matrixSolution):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_zlabel('z axis') 
    for i in range(2):
        if i == 0:
            for i in range(len(matrixOriginal)):
                x1, x2 = matrixOriginal[i][0][0], matrixOriginal[i][1][0]
                y1, y2 = matrixOriginal[i][0][1], matrixOriginal[i][1][1]
                z1, z2 = matrixOriginal[i][0][2], matrixOriginal[i][1][2]
        #        print((x1, y1, z1),(x2, y2, z2))
                plt.plot([x1,x2],[y1,y2],[z1,z2], color = 'blue')
        elif i == 1:
            for i in range(len(matrixSolution)):
                x1, x2 = matrixSolution[i][0][0], matrixSolution[i][1][0]
                y1, y2 = matrixSolution[i][0][1], matrixSolution[i][1][1]
                z1, z2 = matrixSolution[i][0][2], matrixSolution[i][1][2]
        #        print((x1, y1, z1),(x2, y2, z2))
                plt.plot([x1,x2],[y1,y2],[z1,z2], color = 'red')
    plt.show()
    return None

def mainFunc():
    constants = getTemplateAndSolution()
#    print(constants.template, constants.solution, constants.original)
    tuplesOriginal = figureBondedAtoms(constants.template, constants.original)
    tuplesSolution = figureBondedAtoms(constants.template, constants.solution)
    cartPlotMatrixOriginal = gatherCartesianData(tuplesOriginal, constants.original)
    cartPlotMatrixSolution = gatherCartesianData(tuplesSolution, constants.solution)
    plotting(cartPlotMatrixOriginal, cartPlotMatrixSolution)
    return None

mainFunc()
