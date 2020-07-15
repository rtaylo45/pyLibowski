from origenTools.pyTransition import transitionMatrix
from matexp.CRAM import apply
import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
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
nuclideNames = "origenTools/data/baseCaseOrigen.txt"

# Transition matrix object
transition = transitionMatrix(diagDecayFile, diagRxFile, offDiagRxFile)

nuclideDict = {'U': [235],
               'Xe': ['135m', 135],
               'I': [135]
                }
# Set the nuclides
#transition.setProblemNuclides(nuclideDict)
transition.setProblemNuclidesFromFile(nuclideNames)
nuclides = transition.getNuclides()
n0 = np.zeros((len(nuclides), 1))
#n0[:] = 1e20
# initial condition of U 235
for index, nuclideID in enumerate(nuclides):
    isotopeName = transition.getNameFromID(nuclideID)
    if isotopeName == 'U-235':
        n0[index] = 1e10

# Tend
tend = 1000000
steps = 100
dt = tend/steps
t = 0
time = []
solDic = OrderedDict()
# neutron flux
flux = 1.e13
# builds the transition matrix
#matrix = transition.buildTransitionMatrix(flux)

transition.writeLibowskiSpeciesInputFile("speciesInputNames.dat")
transition.writeLibowskiSpeciesReactionFile("speciesInputDecay.dat", decayOnly=True)
transition.writeLibowskiSpeciesReactionFile("speciesInputTrans.dat", transOnly=True)

#transition.writeLibowskiSpeciesInputFile("speciesInputNamesSmall.txt")
#transition.writeLibowskiSpeciesReactionFile("speciesInputDecaySmall.txt", decayOnly=True)
#transition.writeLibowskiSpeciesReactionFile("speciesInputTransSmall.txt", transOnly=True)

"""
plt.spy(matrix)
plt.show()

# Runs the problem
for step in range(steps):
    t = t+dt
    time.append(t)
        
    sol = apply(matrix, dt, n0)
    solDic = unpackSolution(transition, solDic, sol)
    n0 = sol
time = np.asarray(time)
# Loops through the soltuion to plot
for nuclideID in solDic.keys():
    sol = solDic[nuclideID]
    if sol[-1] > 1.:
        isotopeName = transition.getNameFromID(nuclideID)
        plt.plot(time/60./60./24., sol, label=isotopeName)
plt.grid()
plt.legend()
plt.xlabel("Time [day]")
plt.ylabel("Atomic number density")
plt.yscale("log")
plt.show()
"""
