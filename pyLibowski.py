from origenTools.pyTransition import transitionMatrix
from matexp.CRAM import apply
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import sys
from collections import OrderedDict

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
        #oldSol = newDic[nuclide]
        #print(index, nuclide, len(newDic[nuclide]))
        #oldSol.append(sol[index])
        #newDic[nuclide] = oldSol
        newDic[nuclide].append(sol[index])
    return newDic

# Origen files
diagDecayFile = "origenTools/data/diag-dec.txt"
diagRxFile = "origenTools/data/diag-rx.txt"
offDiagRxFile = "origenTools/data/offdiag.txt"
nuclideNames = sys.argv[1]
# file base name
#fbaseName = sys.argv[2]
fbaseName = "moleProblems"

# Transition matrix object
transition = transitionMatrix(diagDecayFile, diagRxFile, offDiagRxFile)


nuclideDict = {'U': [233, 234, 235, 236, 237, 238, 239],
               'Pu':[239],
               'Np':[239],
                }

# Set the nuclides
#transition.setProblemNuclides(nuclideDict)
transition.setProblemNuclidesFromFile(nuclideNames)

# Writes out the files
#transition.writeLibowskiSpeciesInputFile(fbaseName+"SpeciesInputNames.dat")
#transition.writeLibowskiSpeciesReactionFile(fbaseName+"SpeciesInputDecay.dat", decayOnly=True)
#transition.writeLibowskiSpeciesReactionFile(fbaseName+"SpeciesInputTrans.dat", transOnly=True)

transition.writeLibowskiSpeciesInputFile("speciesInputNames.dat")
transition.writeLibowskiSpeciesReactionFile("speciesInputDecay.dat", decayOnly=True)
transition.writeLibowskiSpeciesReactionFile("speciesInputTrans.dat", transOnly=True)

#transMatrix = transition.buildTransitionMatrix(1.e13)
#print(transMatrix)
