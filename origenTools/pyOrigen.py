#print("end", parentNuclideID, daughterNuclideID, type(dfDaughter))
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
        
