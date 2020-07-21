"""
Author: Zack Taylor

transition matrix class that uses the ORIGEN data class to build a trasition
matrix for libowski
"""
import numpy as np
import numpy.linalg as LA
import sys
from .pyOrigen import ORIGENData
from collections import OrderedDict

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

    def setProblemNuclidesFromFile(self, fname):
        """
        Sets the nuclides in the transition matrix from a file

        @param fname    Location of the file
        """
        nuclideDict = self._dataObject.processNuclideList(fname)
        self.setProblemNuclides(nuclideDict)

    def setProblemNuclides(self, nuclideDic, removeDublicateGroups=False):
        """
        Sets the nuclides to be in the transition matrix

        @param nuclideDic   Dictionary that contains the nuclides for the 
                            transition matrix. The 'key' is the element
                            the 'values' is an array of the atomic mass 
                            numbers.
        """
        nuclideIDs = []
        # Loops through the dictionary 
        for atomicNumber, massNumbers in nuclideDic.items():
            if massNumbers:
                for massNumber in massNumbers:
                    IDs = self._dataObject.getNuclideOrigenIDs(atomicNumber, 
                        massNumber, addG1=False)
                    for ID in IDs:
                        nuclideIDs.append(ID)
            else:
                IDs = self._dataObject.getNuclideOrigenIDs(atomicNumber,
                    massNumber)
                for ID in IDs:
                    nuclideIDs.append(ID)

        # Removes dublicate groups from the nuclide list. This part 
        # is kinda sketchy right now. I haven't looked though the file
        # to check if the reaction rates are the same for a nuclide in 
        # different sub groups. 
        # Need to remove this bullshit
        if removeDublicateGroups:
            uniqueNuclides = []
            
            for thisNuclide in nuclideIDs:
                thisGroupID = self._dataObject.getSubGroup(thisNuclide)
                thisVal = int(str(thisNuclide)[1:])
                candidateNuclide = thisNuclide
                dublicateFound = False
                tempDublicateNuclides = []
                for otherNuclide in nuclideIDs:
                    if thisNuclide == otherNuclide:
                        pass
                    else:
                        otherGroupID = self._dataObject.getSubGroup(otherNuclide)
                        otherVal = int(str(otherNuclide)[1:])
                        if thisVal == otherVal:
                            dublicateFound = True
                            tempDublicateNuclides.append(otherNuclide)
                if not dublicateFound:
                    uniqueNuclides.append(candidateNuclide)
                else:
                    if thisGroupID == 3 and dublicateFound:
                        uniqueNuclides.append(candidateNuclide)

            self._nuclides = np.asarray(uniqueNuclides)
        else:
            self._nuclides = np.asarray(nuclideIDs)

    def buildTransitionMatrix(self, flux):
        """
        Builds that transition matrix that is read in my libowski

        @param flux     Neutron flux in 1/cm^2/s

        """
        transMatrix = np.zeros((len(self._nuclides), len(self._nuclides)))
        # builds the matrix index nuclide map
        matrixIndex_nuclide_map = OrderedDict()
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

    def writeLibowskiSpeciesInputFile(self, fname):
        """
        Writes an input file used to set the species. 
        Format example:
            eAName  Mass[lbm/mol]   initial amount[mole]    diffusion coeff[ft^2/s]
            u-235   235.0439        1.0                     0.0

        @param fname    Name of file to write to
        """
        f = open(fname, "w")
        for nuclide in self._nuclides:
            name = self._dataObject.convertIDtoEAmName(nuclide)
            molarMass = str(self._dataObject.getMolarMass(nuclide))
            initCon = str(0.0)
            difCoeff = str(0.0)
            decayConst = str(self._dataObject.getDecayConstant(nuclide))
            string = name + '\t' + molarMass + '\t' + initCon + '\t' + difCoeff + '\n'
            f.write(string)
        f.close()

    def writeLibowskiSpeciesReactionFile(self, fname, decayOnly=False, transOnly=False):
        """
        Writes an input file used to set reaction rates for speices in libowski. In
        libowski the coefficient used for source terms are housed in a std::vector<double> 
        of length numOfSpecies.

        When writting this file it is assumed that the speices ID's in libowsk follow the 
        order that they are from the writeLibowskiSpeciesInputFile function
        """
        f = open(fname, "w")
        # builds the matrix index nuclide map. Not that the index needs to be the same 
        # integer that is returned from the libowski add species function
        matrixIndex_nuclide_map = OrderedDict()
        for index, nuclide in enumerate(self._nuclides):
            matrixIndex_nuclide_map[nuclide] = index
        # Loops though nuclides
        for nuclide in self._nuclides:
            name = self._dataObject.convertIDtoEAmName(nuclide)
            coeffVect = np.zeros((1,len(self._nuclides)))
            # The decay constant 1/s
            diagCoeff = self._dataObject.getDecayConstant(nuclide) 
            thisIndex = matrixIndex_nuclide_map[nuclide]
            if not transOnly:
                coeffVect[0,thisIndex] += -diagCoeff

            # Loops over source terms from neutron induced reactions and
            # decay. These are the off diagonal elements and will be possitive
            for parent in self._dataObject.getReactionParents(nuclide):
                # Loop over parent nuclides
                if parent in self._nuclides:
                    parentIndex = matrixIndex_nuclide_map[parent]
                    # Gets the coefficient. 1/s
                    coeff = self._dataObject.getReactionRate(parent, nuclide, 
                        decayOnly=decayOnly, transOnly=transOnly)
                    # Sets the coefficient
                    coeffVect[0,parentIndex] += coeff
            # the line string to write to the file
            string = str(thisIndex) + '\t' + name + '\t'
            # adds the coefficent to the string
            for index in range(coeffVect.shape[1]):
                coeff = coeffVect[0,index]
                string += np.str(coeff) + '\t'
            string += '\n' 
            f.write(string)
        f.close()

               

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
