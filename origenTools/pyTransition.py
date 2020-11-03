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
    # numpy array of nuclide phase indices 0 = liquid; 1 = gas; 2 = wall
    _phase = None

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

    def _orderIDs(self, nuclideIDs, phases = np.empty(0)):
        """
        Order the nuclide IDs based on their ZAI index
        
        @param  nuclideIDs   Array of nuclide IDs
        @return oIDs         Ordered IDs
        """
        # Get ZAI IDs 
        ZAIs = self._dataObject.getZAI(nuclideIDs)
        if phases.size != 0:
            # Creates array that combines both columns
            array = np.c_[nuclideIDs, ZAIs, phases]
            # sorts the array based on the ZAIs
            array = array[np.argsort(array[:, 1])]
            # Gets the nuclide IDs to set. Ordered by the ZAIs
            oIDs = array[:,0]
            oPhases = array[:,2]
            return oIDs, oPhases
        else:
            # Creates array that combines both columns
            array = np.c_[nuclideIDs, ZAIs]
            # sorts the array based on the ZAIs
            array = array[np.argsort(array[:, 1])]
            # Gets the nuclide IDs to set. Ordered by the ZAIs
            oIDs = array[:,0]
            return oIDs

    def _convertPhaseIndexToPhaseName(self, index):
        """
        Converts the phase index to a string name that will be added to the
        end of a nuclide name
        """
        assert(index <= 2 and index >= 0)
        name = "None"
        if index == 0:
            name = "Liq"
        elif index == 1:
            name = "Gas"
        elif index == 2:
            name = "Wall"
        return name

    def setProblemNuclidesFromExistingArray(self, fname):
        """
        Sets the nuclides from an exsiting file that already has the 
        nuclides in the ORIGEN nuclide format

        @param fname    Location of the file
        """
        array = np.loadtxt(fname, dtype=np.int64)

        nuclideIDs, phases = array[:,0], array[:,1]

        nuclideIDs, phases = self._orderIDs(nuclideIDs, phases)

        self._nuclides = nuclideIDs
        self._phases = phases

    def setProblemNuclidesFromFile(self, fname):
        """
        Sets the nuclides in the transition matrix from a file

        @param fname    Location of the file
        """
        nuclideDict = self._dataObject.processNuclideList(fname)
        self.setProblemNuclides(nuclideDict, False)

    def setProblemNuclides(self, nuclideDic, removeDublicateGroups=True):
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
                        massNumber)
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

            uniqueNuclides = self._orderIDs(np.asarray(uniqueNuclides))
            self._nuclides = self._orderIDs(uniqueNuclides)
            np.savetxt('uniqueNuclides.txt', uniqueNuclides, fmt='%8i')
        else:
            nuclideIDs = self._orderIDs(np.asarray(nuclideIDs))
            self._nuclides = nuclideIDs
            np.savetxt('nuclides.txt', nuclideIDs, fmt='%8i')

    def buildTransitionMatrix(self, flux):
        """
        Builds that transition matrix that is read in my libowski

        @param flux     Neutron flux in 1/cm^2/s

        """
        transMatrix = np.zeros((len(self._nuclides), len(self._nuclides)))
        # builds the matrix index nuclide map
        matrixIndex_nuclide_map = OrderedDict()
        for index, nuclide in enumerate(self._nuclides):
            phaseName = self._convertPhaseIndexToPhaseName(self._phases[index])
            name = self._dataObject.convertIDtoEAmName(nuclide) + phaseName
            matrixIndex_nuclide_map[name] = index
        # Loops though to build the transition matrix
        for index, nuclide in enumerate(self._nuclides):
            thisPhaseName = self._convertPhaseIndexToPhaseName(self._phases[index])
            thisName = self._dataObject.convertIDtoEAmName(nuclide) + thisPhaseName
            thisIndex = matrixIndex_nuclide_map[thisName]

            # The decay constant 1/s
            diagCoeff = self._dataObject.getDecayConstant(nuclide) 
            # Removal rate from neutron induced reactions 1/s
            diagCoeff += self._dataObject.getReactionRemovalRate(nuclide, flux)
            # Sets the coefficient
            transMatrix[thisIndex, thisIndex] = -diagCoeff
           
            # Loops over source terms from neutron induced reactions and
            # decay. These are the off diagonal elements and will be possitive
            for pIndex, parent in enumerate(self._nuclides):
                parentPhaseName = self._convertPhaseIndexToPhaseName(self._phases[pIndex])
                parentName = self._dataObject.convertIDtoEAmName(parent) + parentPhaseName
                # Loop over parent nuclides
                if parent in self._dataObject.getReactionParents(nuclide) and thisPhaseName == parentPhaseName:
                    parentIndex = matrixIndex_nuclide_map[parentName]
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
        for index, nuclide in enumerate(self._nuclides):
            phaseName = self._convertPhaseIndexToPhaseName(self._phases[index])
            name = self._dataObject.convertIDtoEAmName(nuclide) + phaseName
            molarMass = str(self._dataObject.getMolarMass(nuclide))
            initCon = str(0.0)
            difCoeff = str(0.0)
            decayConst = str(self._dataObject.getDecayConstant(nuclide))
            string = name + '\t\t\t' + molarMass + '\t\t\t' + initCon + '\t\t\t' + difCoeff + '\n'
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
            phaseName = self._convertPhaseIndexToPhaseName(self._phases[index])
            name = self._dataObject.convertIDtoEAmName(nuclide) + phaseName
            matrixIndex_nuclide_map[name] = index
        # Loops though nuclides
        for thisIndex, nuclide in enumerate(self._nuclides):
            thisPhaseName = self._convertPhaseIndexToPhaseName(self._phases[thisIndex])
            thisName = self._dataObject.convertIDtoEAmName(nuclide) + thisPhaseName
            thisMM = self._dataObject.getMolarMass(nuclide)
            coeffVect = np.zeros((1,len(self._nuclides)))
            # The decay constant 1/s
            diagCoeff = self._dataObject.getDecayConstant(nuclide) 
            thisIndex = matrixIndex_nuclide_map[thisName]
            if not transOnly:
                coeffVect[0,thisIndex] += -diagCoeff

            # Loops over source terms from neutron induced reactions and
            # decay. These are the off diagonal elements and will be possitive
            for pIndex, parent in enumerate(self._nuclides):
                parentPhaseName = self._convertPhaseIndexToPhaseName(self._phases[pIndex])
                parentName = self._dataObject.convertIDtoEAmName(parent) + parentPhaseName
                parentMM = self._dataObject.getMolarMass(parent)
                # Loop over parent nuclides
                if parent in self._dataObject.getReactionParents(nuclide) and thisPhaseName == parentPhaseName:
                    parentIndex = matrixIndex_nuclide_map[parentName]
                    # Gets the coefficient. 1/s
                    MMratio = thisMM/parentMM
                    coeff = MMratio*self._dataObject.getReactionRate(parent, nuclide, 
                        decayOnly=decayOnly, transOnly=transOnly)
                    # Sets the coefficient
                    coeffVect[0,parentIndex] += coeff
            # the line string to write to the file
            string = str(thisIndex) + '\t' + thisName + '\t'
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
