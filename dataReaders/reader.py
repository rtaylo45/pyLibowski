import numpy as np
from collections import OrderedDict

def read2Ddata(fname):
    returnDict = OrderedDict()
    solverDict = None
    solverName = None
    with open(fname, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith('Solver'):
                candidateSolver = str(line.split()[-1])
                if solverName != candidateSolver:
                    if solverDict:
                        returnDict[solverName] = solverDict
                    solverDict = OrderedDict()
                    solverName = candidateSolver
            elif line.startswith('Total problem time'):
                totalTime = float(line.split()[-1])
            elif line.startswith('xLength'):
                xLength = float(line.split()[-1])
            elif line.startswith('yLength'):
                yLength = float(line.split()[-1])
            elif line.startswith('dt'):
                dt = float(line.split()[-1])
            elif line.startswith('xCells'):
                xCells = float(line.split()[-1])
                numberOfxElements = int(xCells)
            elif line.startswith('yCells'):
                yCells = float(line.split()[-1])
                numberOfyElements = int(yCells)
            elif line.startswith('variables'):
                variables = line.split()[1:len(line.split())]
                elements = numberOfxElements*numberOfyElements
            elif line.startswith('time'):
                time = line.split()[-1]
                solution = np.zeros((elements, len(variables)))
                di = 0
            elif line.startswith('end'):
                returnDict[solverName] = solverDict
                for key in solverDict.keys():
                    data = solverDict[key]
            elif len(line.split()) == 0:
                solverDict[str(time)] = solution
                di = 0
            else:
                splitline = line.split()
                for varIndex in range(len(variables)):
                    solution[di,varIndex] = float(splitline[varIndex])
                di += 1

    return returnDict

def convertTo2Dgrid(numXpoints, numYpoints, data):
    """
    Converts data from libowsk to 2D data grid for seaborn heat map plots

    Converts:           to:
    x1 y1 dat               x1  x2  x3  ... xn
    x1 y2 dat           yn  dat dat dat ... dat
    .  .  .             .
    .  .  .             .
    .  .  .             .
    x1 yn dat           y3  dat dat dat ... dat
    x2 y1 dat           y2  dat dat dat ... dat
    x2 y2 dat           y1  dat dat dat ... dat
    .  .  .
    .  .  .
    .  .  .
    xn y1 dat
    xn y2 dat
    .  .  .
    .  .  .
    .  .  .
    xn yn dat

    x and y are not included.. only the dat array
    """
    array = np.zeros((numYpoints, numXpoints))
    absIndex = 0
    for xIndex in range(numXpoints):
        for yIndex in range(numYpoints-1, -1, -1):
            val = data[absIndex]
            array[yIndex, xIndex] = val
            absIndex += 1

    return array
