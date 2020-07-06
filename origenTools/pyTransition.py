"""
Author: Zack Taylor

transition matrix class that uses the ORIGEN data class to build a trasition
matrix for libowski
"""
import numpy as np
import numpy.linalg as LA
import sys
from pyOrigen import ORIGENData

class transitionMatrix:

    """
    Private data members
    """
    # Origen data object
    _dataObject = None
    # Numpy array of nuclide IDs that will be used to build the transition 
    # matrix
    _nuclides = None

    """
    Class Methods
    """
    def __init__(self, diagDecayFname, diagRxFname, offDiagRxFname):
        """
        Initilizer

        @param diagDecayFname   File path for diagional decay file
        @param diagRxFname      File path for diagional rx file
        @param offDiagRxFname   File path for off diagional rx file
        """
        self._dataObject = ORIGENData(diagDecayFname, diagRxFname, offDiagRxFname)

    def setProblemNuclides(self, nuclideDic, metastable=True):
        """
        Sets the nuclides to be in the transition matrix

        @param nuclideDic   Dictionary that contains the nuclides for the 
                            transition matrix. The 'key' is the element
                            the 'values' is an array of the atomic mass 
                            numbers.
        @param metastable   Set to true if the user wants to include 
                            metastable nuclides 
        """
        nuclideIDs = []
        # Loops through the dictionary 
        for atomicNumber, massNumbers in nuclideDic.items():
            if massNumbers:
                for massNumber in massNumbers:
                    IDs = self._dataObject.getNuclideOrigenIDs(atomicNumber, 
                        massNumber, addG1=False)
                    for ID in IDs:
                        if self._dataObject.isMetastable(ID) and not metastable:
                            pass
                        else:
                            nuclideIDs.append(ID)
            else:
                IDs = self._dataObject.getNuclideOrigenIDs(atomicNumber,
                    massNumber, addG1=False)
                for ID in IDs:
                    nuclideIDs.append(ID)
        self._nuclides = np.asarray(nuclideIDs)

    def buildTransitionMatrix(self, flux):
        """
        Builds that transition matrix that is read in my libowski

        @param flux     Neutron flux in 1/cm^2/s

        """
        transMatrix = np.zeros((len(self._nuclides), len(self._nuclides)))
        # builds the matrix index nuclide map
        matrixIndex_nuclide_map = {}
        for index, nuclide in enumerate(self._nuclides):
            matrixIndex_nuclide_map[nuclide] = index

        # Loops though to build the transition matrix
        for nuclide in self._nuclides:
            thisIndex = matrixIndex_nuclide_map[nuclide]

            # The decay constant 1/s
            diagCoeff = self._dataObject.getDecayConstant(nuclide) 
            # Removal rate from neutron induced reactions 1/s
            diagCoeff += self._dataObject.getReactionRemovalRate(nuclide, flux)
            # Sets the coefficient
            transMatrix[thisIndex, thisIndex] = -diagCoeff
           
            # Loops over source terms from neutron induced reactions and
            # decay. These are the off diagonal elements and will be possitive
            for parent in self._dataObject.getReactionParents(nuclide):
                # Loop over parent nuclides
                if parent in self._nuclides:
                    parentIndex = matrixIndex_nuclide_map[parent]
                    # Gets the coefficient. 1/s
                    coeff = self._dataObject.getReactionRate(parent, nuclide, flux)
                    # Sets the coefficient
                    transMatrix[thisIndex, parentIndex] += coeff

        return transMatrix

    def getNuclides(self):
        """
        Retunes a list of the nuclides in nuclideID form that are in the 
        transition matrix. This is in the same order as the solution vector.
        """
        return self._nuclides

    def write(self, fname):
        """
        Writes the transition matrix to an output file.

        """
        pass

    def getNameFromID(self, nuclideID):
        """
        Returns the nuclideID name of an isotope in EAm form:
            E - Element i.e. Xe
            A - Mass number 
            m - Metastable indecator
        The returned name for the ID of Xenon 135 would be 
        Xe-135

        @param nuclideID    Origen ID of the nuclide
        """
        return self._dataObject.convertIDtoEAmName(nuclideID)


if __name__=="__main__":
    #from pyLibowski import CRAM
    import matplotlib.pyplot as plt

    def unpackSolution(transObj, solDic, sol):
        nuclides = transObj.getNuclides()
        newDic = solDic
        # if the solDic is empty (it will be when the first solution
        # is upbacked) it builds the solution Dict
        if not newDic:
            for nuclide in nuclides:
                newDic[nuclide] = []
        # loops over nuclides to update solution
        for index, nuclide in enumerate(nuclides):
            oldSol = newDic[nuclide]
            oldSol.append(sol[index])
            newDic[nuclide] = oldSol
        return newDic

    diagDecayFile = "data/diag-dec.txt"
    diagRxFile = "data/diag-rx.txt"
    offDiagRxFile = "data/offdiag.txt"
    flux = 1e13

    transition = transitionMatrix(diagDecayFile, diagRxFile, offDiagRxFile)

    nuclideDic = {'Xe': ['135m'],
                  'I' : ['135'],
                  #'Sb': [135],
                  #'Nd': [149],
                  #'Te': [135],
                  'U' : ['238']
                  }
    transition.setProblemNuclides(nuclideDic, metastable=True)
    nuclides = transition.getNuclides()
    print(nuclides)
    n0 = np.zeros((len(nuclides), 1))
    n0[-1] = 1e10
    tend = 500*100
    steps = 10000
    dt = tend/steps
    t = 0
    time = []
    solDic = {}
    matrix = transition.buildTransitionMatrix(0)
    print (matrix)
   
    for step in range(steps):
        t = t+dt
        time.append(t)
        if t < 1000:
            flux = 1e13
        else:
            flux = 0.0
            
        matrix = transition.buildTransitionMatrix(flux)
        sol = CRAM.apply(matrix, dt, n0)
        solDic = unpackSolution(transition, solDic, sol)
        n0 = sol
    """
    flux = 1e13
    matrix = transition.buildTransitionMatrix(flux)
    for step in range(steps):
        t = t+dt
        time.append(t)
        sol = CRAM.apply(matrix, t, n0)
        solDic = unpackSolution(transition, solDic, sol)
    """
  
    for nuclideID in solDic.keys():
        sol = solDic[nuclideID]
        isotopeName = transition.getNameFromID(nuclideID)
        plt.plot(time, sol, label=isotopeName)
    plt.grid()
    plt.legend()
    plt.xlabel("Time [s]")
    plt.ylabel("Atomic number density")
    plt.yscale("log")
    plt.show()
