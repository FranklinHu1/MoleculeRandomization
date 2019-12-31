#generates histograms
from errorMinSubRoutine import run
import matplotlib.pyplot as plt
import numpy as np

'''
Module used to show variation in interatomic distances. Used during the starting
stages of the project to show that there was a uniform sampling of bond lengths.
No longer necessary for the project to proceed.
'''

def generateData(times):
    x = []
    interestedString = 'CC'
    reverseString = 'CC'
    iterations = 0
    for i in range(times):
        tempMatrix = run()[0]
        for row in range(len(tempMatrix)):
            checkString = ''
            for c in tempMatrix[row,0]:
                if c.isalpha():
                    checkString += c
            if checkString == interestedString or checkString == reverseString:
                x.append(tempMatrix[row,1])
        iterations += 1
        print(iterations)
    datums = np.array(x, dtype = float)
    return datums
    
def generateHistogram(times):
    data = generateData(times)
    plt.hist(data, 50)
    plt.xlabel('C-C distances in Angstroms')
    plt.ylabel('Frequency')
    plt.title('Histogram of C-C distances, data generated ' + str(times) + ' times')
    
generateHistogram(10000) #The molecule right now is C6N1H11, angle is +- 10 degrees
 
#x = [1,3,4,6,7,3,3,2,4,6,7,5,3,4,5,1,3,4,6,7,8,9,5,3,2,4,5,6,1,3,4,4]
#plt.hist(x, [1,2,3,4,5,6,7,8])    
#        
#
#
#x = [1.7009406230355648 ,1.367162790770352, 1.2561705927929918, 1.59463076063187,
# 1.2291658082070906]
#
#plt.hist(x, [0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8])
