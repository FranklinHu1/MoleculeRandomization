#Specifically works for qm7
import scipy.io as sio
from mendeleev import element
import numpy as np
from copy import deepcopy
import xlwt
'''
Early code for extracing molecular data from the qm7 dataset. This was early on
when the project was in its nascent phases, and we were doing a proof-of-concept
type thing. No longer necessary for the project.
'''

fullPath = 'C:\\Users\Frank\Box Sync\College\Research and Additional Folders\MoleculeRandomizationUpdated\Molecules\\'
#You will want to change the path name for your specific directories
moleculeName = 'C6N1H11'
#Change this string for the desired molecule

def readMatLabData(dataset, moleculeRow):   
    mat_contents = sio.loadmat(dataset)
        
    masterAtomMatrix = mat_contents['Z']
    masterCartesianMatrix = mat_contents['R']
    
    masterAtomMatrix = np.array(masterAtomMatrix, dtype = object)
    masterCartesianMatrix = np.array(masterCartesianMatrix, dtype = object)
#changes the heading rows in masterAtomMatrix to the actual atoms
    stopPoint1, stopPoint2 = None, None
    for i in range(len(masterAtomMatrix[0])):
        if masterAtomMatrix[moleculeRow, i] == float(0):
            stopPoint1 = i
            break
    nonZeros = masterAtomMatrix[moleculeRow,:stopPoint1]
    for i in range(len(nonZeros)):
        if nonZeros[i] == float(1):
            stopPoint2 = i
            break
    nonZerosAndNonOnes = nonZeros[:stopPoint2]
    nonZeros = sorted(nonZerosAndNonOnes) + [1.0]*((len(nonZeros))-stopPoint2)
    temp = nonZeros + [0.0]*((len(masterAtomMatrix[0]))-stopPoint1)
    masterAtomMatrix[moleculeRow] = np.array(temp, dtype = object)
    for col in range(len(masterAtomMatrix[0])):
        atomicNumber = int(masterAtomMatrix[moleculeRow,col])
        if atomicNumber != 0:
            masterAtomMatrix[moleculeRow,col] = str(element(atomicNumber).symbol) + str(col+1)
        else:
            pass
    return [masterAtomMatrix, masterCartesianMatrix]


def matrixEditing(matrix):
    temporary = deepcopy(matrix)
    editedMat = []
    headerList, xList, yList, zList = [],[],[],[]
    for i in range(len(temporary[0,0])):
        if temporary[0,0][i] != float(0) and type(temporary[0,0][i]) == str:
            headerList.append(temporary[0,0][i])
        if temporary[1,0][i] != float(0) and type(temporary[1,0][i]) == float:
            xList.append(temporary[1,0][i])
        if temporary[2,0][i] != float(0) and type(temporary[2,0][i]) == float:
            yList.append(temporary[2,0][i])
        if temporary[3,0][i] != float(0) and type(temporary[3,0][i]) == float:
            zList.append(temporary[3,0][i])
    counter = 1
    for k in range(len(headerList)):
        if k != (len(headerList)-1):
            if headerList[k][0] == headerList[k+1][0]:
                headerList[k] = headerList[k][0] + str(counter)
                counter += 1
            elif headerList[k][0] != headerList[k+1][0]:
                headerList[k] = headerList[k][0] + str(counter)
                counter = 1
        else:
            headerList[k] = headerList[k][0] + str(counter)
    editedMat.append(headerList)
    editedMat.append(xList)
    editedMat.append(yList)
    editedMat.append(zList)
    editedMat = np.array(editedMat, dtype = object)
    return editedMat


def generateMolecule(moleculeRow, masterAtomMatrix, masterCartesianMatrix):
    
    atoms = masterAtomMatrix[moleculeRow,:]
    xCoords = masterCartesianMatrix[moleculeRow,:,0]
    yCoords = masterCartesianMatrix[moleculeRow,:,1]
    zCoords = masterCartesianMatrix[moleculeRow,:,2]
    
    finalMat = np.empty((4,1), dtype = object)
    finalMat[0,0] = atoms
    finalMat[1,0] = xCoords
    finalMat[2,0] = yCoords
    finalMat[3,0] = zCoords
    
    finalMat = matrixEditing(finalMat)
    
    return finalMat

def generateWorkbookName(matrix):
    copyMat = deepcopy(matrix)
    fileName = ''
    atomCountDict = dict()
    for col in range(len(copyMat[0])):
        atomHeader = copyMat[0][col]
        atomAndNumList = atomHeader.split(atomHeader[1]) #splits to get a list of atom and number
        if atomAndNumList[0] in atomCountDict:
            atomCountDict[atomAndNumList[0]] += 1
        else:
            atomCountDict[atomAndNumList[0]] = 1
    for key in atomCountDict:
        fileName += key + str(atomCountDict[key])
    return fileName

def generateExcelFile(matrix, fileName):
    wb = xlwt.Workbook()
    ws = wb.add_sheet('Molecule')
    for row in range(len(matrix)):
        if matrix[row] != []:
            for col in range(len(matrix[0])):
                ws.write(row, col, matrix[row][col])
        else:
            pass
    wb.save(fullPath + fileName + '.xls')
    return None
    
def mainFunc(dataSet, moleculeRow):
    rawMatrixOutput = readMatLabData(dataSet, moleculeRow)
    moleculeMat = generateMolecule(moleculeRow, rawMatrixOutput[0], rawMatrixOutput[1])
    fileNames = generateWorkbookName(moleculeMat)
    generateExcelFile(moleculeMat, fileNames)
    return None

mainFunc('qm7.mat',1)
    










        
        
        
        

        





