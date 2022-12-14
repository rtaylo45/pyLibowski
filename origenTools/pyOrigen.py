"""
Author: Zack Taylor
Origen data class that does post processing of SCALE transition matrix data
obtained using obiwan. To get these matrices from origen:
obiwan view -type=coefft path/to/origen.f33 >offdiag.txt
obiwan view -type=loxs path/to/origen.f33 >diag-rx.txt
obiwan view -type=nucl path/to/origen.f33 >diag-dec.txt
Using these data files, the ORIGENData class builds an input file for libowski
to use for setting up depletion problems.
"""
from collections.abc import Iterable
import numpy as np
import pandas as pd
from periodictable import elements
import sys
from collections import OrderedDict

class ORIGENData:

    """
    Private data members
    """
    # Pandas array that holds the diagonal reaction data for
    # the transition matirx. Order of the Data.
    # nuclide, mass[g/mol], abondence[atom%], decay[1/s], reaction[barns]
    _diagionalPandasData = None
    # Pandas array that holds the off diagonal reaction data for
    # the transition matrix
    _offDiagionalPandasData = None

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
        # Processes the diagional data files and builds the
        # diagionalPandasData variable
        self._diagionalPandasData = self._processDiagionalData(diagDecayFname,
            diagRxFname)
        # Processes the off diagional data file
        self._offDiagionalPandasData = self._processOffDiagionalData(offDiagRxFname)

    def _processDiagionalData(self, diagDecayFname, diagRxFname):
        """
        Process the diagional data files
        @param diagDecayFname   File path for diagional decay file
        @param diagRxFname      File path for diagional rx file
        @return     Pandas data frame of diagional data
        """
        # reads the diagional decay data file
        diagDataFrame = pd.read_csv(diagDecayFname, sep='\s+|\t+', engine='python')
        # reads the diagional removal rx data file
        diagRxDataFrame = pd.read_csv(diagRxFname, sep='\s+|\t+', engine='python')
        # drops the nuclide data column from the rx data frame because it is
        # redundent for the final data frame
        diagRxDataFrame = diagRxDataFrame.drop(columns=['nuclide'])
        # combins the two data frames
        diagDataFrame = pd.concat([diagDataFrame, diagRxDataFrame], axis=1)
        # renames the last column
        diagDataFrame = diagDataFrame.rename(columns={'0.0000e+00': 'rx[barn]'})

        nuclideIDs = diagDataFrame.get('nuclide').to_numpy()
        # Builds the elements data frame
        elements = [self.getElement(nuclideID) for nuclideID in nuclideIDs]
        elements = pd.DataFrame(data=elements, index=None, columns=['atomic number'])
        # Builds the mass number data frame
        massNumbers = [self.getMassNumber(nuclideID) for nuclideID in nuclideIDs]
        massNumbers = pd.DataFrame(data=massNumbers, index=None, columns=['mass number'])
        # combins the data frames
        diagDataFrame = pd.concat([diagDataFrame, elements, massNumbers], axis=1)

        return diagDataFrame

    def _processOffDiagionalData(self, offDiagRxFname):
        """
        Process the off diagional data file
        @param offDiagRxFname   File path for off diagional rx file
        @return     Pandas data frame of off diagional reaction rates
        """
        # reads in the off diagional data file
        offDiagRxDataFrame = pd.read_csv(offDiagRxFname, sep='<--\[|\]--|\s+|\t+',
            engine='python')
        # renames the last column
        offDiagRxDataFrame = offDiagRxDataFrame.rename(columns={'0.0000e+00': 'rx[barn]'})

        return offDiagRxDataFrame

    def processNuclideList(self, fname):
        """
        Processes a nuclide input list and returns a python dictionary that 
        is the nuclide input list 
        @param fname    File name in EAm form
        """
        nuclideDict = OrderedDict()
        # Opens the file and reads the lines
        with open(fname, 'r') as f:
            lines = f.readlines()
            # Loops through lines and splits them
            for line in lines:
                splitline = line.split('-')
                ele = splitline[0]
                loc = splitline[1]
                mass = loc.split()[0]
                if ele not in nuclideDict.keys():
                    nuclideDict[ele] = [mass]
                else:
                    if mass not in nuclideDict[ele]:
                        nuclideDict[ele].append(mass)
        return nuclideDict

    def getZAI(self, nuclideIDs):
        """
        Converts the IZZZAAA Origen index to ZAI index defined by

        ZAI = 10,000*Z + 10*A + I

        Can accept single ID or a np array maybe list?

        @param nuclideID    Nuclide ID based on Origen indexing
        @return ZAI         Array or integer of the converted
                            IZZZAAA ID to ZAI
        """
        if isinstance(nuclideIDs, Iterable):
            ZAIs = []
            for nuclideID in nuclideIDs:
                Z = self.getElement(nuclideID)
                A = self.getMassNumber(nuclideID)
                I = self.getIsomericState(nuclideID)
                ZAI = 10000*Z + 10*A + I
                ZAIs.append(ZAI)
            return np.asarray(ZAIs, dtype=np.int64)
        else:
            Z = self.getElement(nuclideIDs)
            A = self.getMassNumber(nuclideIDs)
            I = self.getIsomericState(nuclideIDs)
            ZAI = 10000*Z + 10*A + I
            return np.int(ZAI)

    def getElement(self, nuclideID):
        """
        Gets the atomic number from the nuclide ID
        @param nuclideID    Nuclide ID based on Origen indexing
        @return atomicID   Atomic number that identifies an element
        """
        # convert to string to pluck out the atomic number easier
        nuclideIDstr = np.str(nuclideID)
        assert(len(nuclideIDstr)==8)
        # get atomic number
        atomicIDstr = nuclideIDstr[2:5]

        return np.int(atomicIDstr)

    def getMassNumber(self, nuclideID):
        """
        Gets the atomic mass number from the nuclide ID
        @param nuclideID    Nuclide ID based on Origen indexing
        @return massID      Atomic mass number
        """
        # convert to string to pluck out the atomic number easier
        nuclideIDstr = np.str(nuclideID)
        assert(len(nuclideIDstr)==8)
        # gets the mass number
        massID = nuclideIDstr[5:]

        return np.int(massID)

    def getIsomericState(self, nuclideID):
        """
        Returns the metastable ID

        @param nuclideID    Nuclide ID based on Origen indexing
        @return             Isomeric state ID
        """
        # convert to string to pluck out the atomic number easier
        nuclideIDstr = str(nuclideID)
        assert(len(nuclideIDstr)==8)
        # gets the isometric state
        metastableID = int(nuclideIDstr[1])

        return np.int(metastableID)

    def isMetastable(self, nuclideID):
        """
        Returns true if the nuclide ID is for a metastable nuclide
        @param nuclideID    Nuclide ID based on Origen indexing
        @return logic       Logical for the metastable nuclide
        """
        metastableID = self.getIsomericState(nuclideID)

        if (metastableID == 0):
            return False
        else:
            return True

    def getSubGroup(self, nuclideID):
        """
        Gets the subgroup number from the nuclide ID
        @param nuclideID    Nuclide ID based on Origen indexing
        @return groupID     The subgroup ID based on Origen logic
        """
        # convert to string to pluck out the atomic number easier
        nuclideIDstr = str(nuclideID)
        assert(len(nuclideIDstr)==8)
        # gets the subgroup number
        groupID = nuclideIDstr[0]

        return int(groupID)

    def getDecayConstant(self, nuclideID):
        """
        Gets the decay constant
        @param nuclideID    Nuclide ID based on Origen indexing
        @return lambda_     Decay constant 1/s
        """
        assert(len(str(nuclideID))==8)
        data = self._diagionalPandasData.set_index('nuclide')
        lambda_ = data.at[nuclideID,'decay[1/s]']
        return lambda_

    def getMolarMass(self, nuclideID):
        """
        Gets the molar mass
        @param nuclideID    Nuclide ID based on Origen indexing
        @return mm          Molar mass [g/mol]
        """
        assert(len(str(nuclideID))==8)
        data = self._diagionalPandasData.set_index('nuclide')
        mm = data.at[nuclideID,'mass[g/mol]']
        return mm

    def getReactionRemovalRate(self, nuclideID, flux=1.0):
        """
        Gets the reaction removal rate in barns
        @param nuclideID    Nuclide ID based on Origen indexing
        @param flux         Neutron flux in 1/cm^2/s
        @return rxRate      Reaction removal rate for nuclide
                            in 1/s
        """
        assert(len(str(nuclideID))==8)
        data = self._diagionalPandasData.set_index('nuclide')
        rxRate = data.at[nuclideID,'rx[barn]']
        return rxRate*1.e-24*flux

    def getReactionRate(self, parentNuclideID, daughterNuclideID, flux=1.0, decayOnly=False,
            transOnly=False):
        """
        Gets the reaction rate in barns from parent nuclide
        to daughter nuclide
        @param parentNuclideID      Parent nuclide ID
        @param daughterNuclideID    Daughter nuclide ID
        @param flux                 Neutron flux in 1/cm^2/s
        @param decayOnly            Logical set to true if the user only wants decay reactions
        #param transmuationOnly     Logical set to tru if the user only wants transmutation
                                    reactions
        @return rxRate              Reaction rate from parent nuclide to
                                    daughter nuclide in 1/s
        """
        dfDaughter = self._offDiagionalPandasData.set_index('daughter').loc[[
            daughterNuclideID]]
        dfReaction = dfDaughter.set_index('parent').loc[parentNuclideID]
        # the parent daughter relation only has multiple reaction 
        # mechanisims
        if isinstance(dfReaction['tid'], Iterable):
            rxnRate = 0.0
            for rxnIndex, reactionID in enumerate(dfReaction['tid']):
                rxn = dfReaction.iloc[rxnIndex]
                # If the reaction is not decay.
                if reactionID != -1:
                    if decayOnly:
                        rxnRate += 0.0
                    else:
                        rxnRate += rxn['rx[barn]']*1.e-24*flux
                else:
                    if transOnly:
                        rxnRate += 0.0
                    else:
                        # The reaction is decay
                        rxnRate += rxn['rx[barn]']
            return rxnRate
        else:
            rxnRate = 0.0
            # If the reaction is not decay.
            # The parent daughter relation only has one reaction 
            if dfReaction['tid'] != -1:
                if decayOnly:
                    rxnRate = 0.0
                else:
                    rxnRate = dfReaction['rx[barn]']*1.e-24*flux
            else:
                # The reaction is decay
                if transOnly:
                    rxnRate = 0.0
                else:
                    rxnRate = dfReaction['rx[barn]']
            return rxnRate

    def getReactionDaughters(self, parentNuclideID, tol=0.0):
        """
        Gets an array of all the daughter nuclides that are generated
        from a nuclear reaction of parent nuclide.
        @param parentNuclideID  Parent nuclide ID
        @param tol              Tolarence of the reaction rate. Only returns
                                daughters that have reaction rate above the
                                tol. If tol is 0.0, returns all daughters.
        @return numpyDaughters  Numpy array of all daughter nuclide ID's
        """
        # gets dataFrame of parent nuclide
        dfParent = self._offDiagionalPandasData.set_index('parent').loc[
            parentNuclideID]
        try:
            numpyDaughters = dfParent['daughter'].to_numpy()
        except:
            numpyDaughters = np.asarray([int(dfParent['daughter'])])
        return numpyDaughters

    def getReactionParents(self, daughterNuclideID, tol=0.0):
        """
        Gets an array of all the parent nuclides that are responsible for
        generation of the daughter nuclide
        @param daughterNuclideID    Daughter nuclide ID
        @param tol                  Tolarence of the reaction rate. Only returns
                                    parents that have reaction rate above the
                                    tol. If tol is 0.0, returns all parents.
        @return numpyParents        Numpy array of all daughter nuclide ID's
        """
        dfDaughter = self._offDiagionalPandasData.set_index('daughter').loc[
            daughterNuclideID]
        try:
            numpyParents = dfDaughter['parent'].to_numpy()
        except:
            numpyParents = np.asarray([int(dfDaughter['parent'])])
        return numpyParents

    def convertIDtoEAmName(self, nuclideID):
        """
        Returns the nuclideID name of an isotope in EAm form:
            E - Element i.e. Xe
            A - Mass number 
            m - Metastable indecator
        The returned name for the ID of Xenon 135 would be 
        Xe-135
        @param nuclideID    Origen ID of the nuclide
        """
        ele = str(elements._element[self.getElement(nuclideID)])
        mass = str(self.getMassNumber(nuclideID))
        name = ele+'-'+mass
        if self.isMetastable(nuclideID):
            name = name +'m'
        return name

    def getNuclideOrigenIDs(self, nuclideElementInput, massNumber=None, addG1=True,
        addG2=True, addG3=True):
        """
        Returns a numpy array of the origen nuclide IDs that match the
        nuclide element input, either the element atomic symbol or the
        atomic number.
        @param nuclideElementInput      Element input, either a string
                                        of the element's symbol or the
                                        atomic number
        @param massNumber               Optional argument that allows
                                        the user to specifiy the mass
                                        number of the element as well. 
                                        This includes the m tag at the 
                                        end to indecate the metastable 
                                        isotope.
        @param addG1                    To include to exclude sub
                                        group one nuclides
        @param addG2                    To include to exclude sub
                                        group two nuclides
        @param addG3                    To include to exclude sub
                                        group three nuclides
        @return numpyOrigenNuclideIDs   An array of the Origen nuclide
                                        ID's that match the element
        """
        # Flag to indecate if the nuclide is a metastable one
        isMetastable = False
        # checks to see if input is the elements symbol or not
        if isinstance(nuclideElementInput, str):
            ele = elements.symbol(nuclideElementInput).number
        else:
            ele = nuclideElementInput
        IDs = []
        # Tried to convert the string mass number to an integer
        try:
            massNumber = int(massNumber)
        # If it can't then than means that the mass number as the 
        # m flag at the end indecating that it is metastable
        except:
            isMetastable = True
            massNumber = int(massNumber[:-1])
        # loops through nuclide ID's
        for nuclideID in self._diagionalPandasData.get('nuclide'):
            # if the nuclide ID's atomic number matches the input
            if self.getElement(nuclideID) == ele:
                # if the mass number is entered only add nuclide ID's that
                # match the mass number
                if massNumber:
                    if massNumber == self.getMassNumber(nuclideID) and isMetastable == self.isMetastable(nuclideID):
                        IDs.append(nuclideID)
                else:
                    IDs.append(nuclideID)

        # removes group 1 nuclides if needed
        if not addG1:
            for nuclide in reversed(IDs):
                if self.getSubGroup(nuclide) == 1:
                    IDs.remove(nuclide)
        # removes group 2 nuclides if needed
        if not addG2:
            for nuclide in reversed(IDs):
                if self.getSubGroup(nuclide) == 2:
                    IDs.remove(nuclide)
        # removes group 3 nuclides if needed
        if not addG3:
            for nuclide in reversed(IDs):
                if self.getSubGroup(nuclide) == 3:
                    IDs.remove(nuclide)


        numpyOrigenNuclideIDs = np.asarray(IDs)
        return numpyOrigenNuclideIDs

    def keepActinide(self, nuclideInput, atomicNumber, massNumber):
        """
        Removes all actinides that are not attributed to the specific
        atomicNumber or massNumber that is entered. Basicialy removes
        all actinides but one
        @param nuclideInput         Array of nuclideIDs
        @param atomicNumber         Atomic number of actinide to keep
        @param massNumber           Mass number of actinide to keep
        @return numpyNuclideOutput  Array of clean nuclide ID's
        """
        nuclideInputList = list(nuclideInput)
        for nuclide in reversed(nuclideInputList):
            if self.getSubGroup(nuclide) == 2:
                nuclideMass = self.getMassNumber(nuclide)
                nuclideEle = self. getElement(nuclide)
                if (nuclideMass == massNumber and nuclideEle == atomicNumber):
                    pass
                else:
                    nuclideInputList.remove(nuclide)

        numpyNuclideOutput = np.asarray(nuclideInputList)
        return numpyNuclideOutput
     
