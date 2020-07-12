from origenTools.pyTransition import transitionMatrix
from origenTools.pyOrigen import ORIGENData

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


# Origen files
diagDecayFile = "origenTools/data/diag-dec.txt"
diagRxFile = "origenTools/data/diag-rx.txt"
offDiagRxFile = "origenTools/data/offdiag.txt"
nuclideNames = "origenTools/data/baseCaseOrigen.txt"

# Transition matrix object
transition = transitionMatrix(diagDecayFile, diagRxFile, offDiagRxFile)
