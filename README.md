# MoleculeRandomization
Code for randomizing and optimizing molecule geometries based on qm7 dataset, using scipy.optimize.least_squares with the Trust Region Reflective algorithm.


# Important Modules and Folders:
## qm7.mat
The raw dataset containing the atomic numbers and cartesian coordinates of each atom in each of 7165 molecules.
## matlabProcessing.py
Pulls molecule data from the qm7.mat dataset, giving the cartesian coordinates (x,y,z) of each atom in the molecule in the form of an excel file. Doing so requires knowing the row index of the molecule in the qm7 file. This module also has the fullPath variable, which directs where molecule cartesian data excel files are saved and pulled from.
## histoGen.py
Generates a histogram of all the inter-atomic distances (e.g. all C-C distances) from a large number of randomizations and optimizations of the molecular geometry. The number of times this data is generated can be set by the user.
## targetOutputBase.py
Processes the molecule data from the molecule excel files, and generates the beginning bonds and angles that are found in the original molecule. The output from this module is set as the starting point for the optimization process.
## errorMinSubRoutine.py
Module that randomizes bond lengths and angles of the molecule and generates a target molecule. Also runs the scipy.optimize.least_squares algorithm for optimizing the structure through minimizing residuals
## targetOutputIterate.py 
Modified version of targetOutputBase.py that operates within the run function of errorMinSubRoutine.py
## moleculePlotting.py
Module for visualizing differences between original molecule and randomized molecule. Still a work in progress.
## Molecules (Folder)
Contains all of the molecule cartesian data excel files. The files are named in accordance to the formula of the molecule (e.g. C2H6 for ethane).
## Histograms (Folder)
Contains the histograms for C-C distances based on different number of runs, different tolerances, different angle and bond variation ranges. Files in this folder are named by the following convention: 
(bond length variation range (angstroms), angle variation range (radians), tolerance for least_squares, number of runs, identifier)
For example, a file name of (0.3, 0.1745, 1e-3, 10k, 2) means that the bonds were varied by +/- 0.3 angstroms, angles were varied by +/- 0.1745 radians, the tolerance was 1e-3 for the least_squares algorithm, the data resulted from 10000 runs, and it is the second histogram of this kind.


# Things to do before getting started
1. Update the fullPath variable in matlabProcessing.py so that it points to the Molecules folder correctly. The other modules import the fullPath variable from this module.
2. Make sure that you are randomizing the correct molecule by checking the molecule name in both targetOutputBase.py and targetOutputIterate.py. Ensure that both modules have the same molecule. 
3. Make sure that the variation ranges for bond lengths and angles is satisfactory. This can be found on lines 64 and 65 of the errorMinSubRoutine.py module, where line 64 deals with angle variation and line 65 deals with bond length variation.
4. Check that the tolerances are satisfactory. This is on line 86 of errorMinSubRoutine.py:                                       ```psol = least_squares(residualCalc,x0, ftol = 1e-3, xtol = 1e-3, gtol = 1e-3) ```Note: decreasing the tolerances (e.g. 1e-3 --> 1e-8) will greatly increase run time of the algorithm, but also improve accuracy of optimization
5. You can set the number of algorithm runs you want in the histoGen.py module, on line 32. The input parameter for the function generateHistogram is the number of times.
6. You can also adjust which inter-atom distances you want. This is on line 9 of histoGen.py, governed by the variable interestedString. The string 'CC' refers to the distance between carbon atoms only. Thus, 'NC' would refer to the distance between nitrogen and carbon atoms. Be sure to also enter the atom pair in reverse order on line 10 for the variable reverseString. 
7. The number of bins for the histogram is set under the function generateHistogram in the histoGen.py module, on line 28. The title for the histogram must also be manually adjusted on lines 29-31 under that function.


# Important things to keep in mind
1. The block comment sections in the errorMinSubRoutine.py module is used for generating a solution matrix, which is composed of the new cartesian coordinates for each atom in the molecule. Uncomment these if you want a solution for an optimization, and add the variable 'solution' to the return statement of the run() function.
2. The console will print out the number of trials completed, so you can track the run progress.
