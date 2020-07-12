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
    print()
    print()
    print()
    return newDic

# Origen files
diagDecayFile = "origenTools/data/diag-dec.txt"
diagRxFile = "origenTools/data/diag-rx.txt"
offDiagRxFile = "origenTools/data/offdiag.txt"
nuclideNames = "origenTools/data/baseCaseOrigen.txt"

# Transition matrix object
transition = transitionMatrix(diagDecayFile, diagRxFile, offDiagRxFile)

# Set the nuclides
#transition.setProblemNuclidesFromFile(nuclideNames)
transition.setProblemNuclidesFromFile(nuclideNames)
nuclides = transition.getNuclides()
print('transition nuclide size', len(nuclides))
n0 = np.zeros((len(nuclides), 1))
# initial condition of U 238
n0[147] = 1e10
# Tend
tend = 500*1000
steps = 100
dt = tend/steps
t = 0
time = []
solDic = OrderedDict()
# neutron flux
flux = 1.e13
# builds the transition matrix
matrix = transition.buildTransitionMatrix(flux)
plt.spy(matrix)
plt.show()

# Runs the problem
for step in range(steps):
    t = t+dt
    time.append(t)
        
    sol = apply(matrix, dt, n0)
    print(sol)
    solDic = unpackSolution(transition, solDic, sol)
    n0 = sol

# Loops through the soltuion to plot
for nuclideID in solDic.keys():
    sol = solDic[nuclideID]
    isotopeName = transition.getNameFromID(nuclideID)
    plt.plot(time, sol, label=isotopeName)
plt.grid()
#plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Atomic number density")
plt.yscale("log")
plt.show()
