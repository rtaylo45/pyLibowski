"""
Author: Zack Taylor

Class that lets me plot things in libowski
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from dataReaders.reader import convertTo2Dgrid
import numpy as np

def plotSingleCellDepletion(fname):
    """
    Plots depletion in a single cell

    @param fname    File name to of the libowski output data
    """
    time, data = pyPlotter.readSingleCellDepletionFile(fname)
    fig = plt.figure()

    ax = plt.subplot(111)
    for specName in data.keys():
        sol = data[specName]
        if sol[-1] > 1.:
            ax.plot(time/60./60./24., sol, label=specName)
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.3, box.width, box.height * 0.8])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), 
            fancybox=True, shadow=True, ncol=15, prop={'size': 8})
    plt.grid()
    plt.xlabel("Time [day]")
    plt.ylabel("Atomic number density")
    plt.yscale("log")
    plt.show()
    
def readSingleCellDepletionFile(fname):
    """
    Reads a the output file for a single cell depletion case

    @param fname    File name of the libowski output data
    """
    time = []
    dataDict = {}
    with open(fname, 'r') as f:
        lines = f.readlines()
        for line in lines:
            splitline = line.split()
            if line.startswith("Time step"):
                time.append(float(splitline[2]))
            elif len(splitline) == 0:
                pass
            else:
                name = splitline[0]
                if name not in dataDict.keys():
                    dataDict[name] = []
                    dataDict[name].append(float(splitline[1]))
                else:
                    dataDict[name].append(float(splitline[1]))
    return np.asarray(time), dataDict


def make2DGIF(data, maxVals):
    fig = plt.figure()
    #ax = plt.axes(ylim=(0, maxVals[0]))
    ax = plt.axes()
    line, = ax.plot([], [])

    ax.set_xlabel('Concentration')
    ax.set_ylabel('y (cm)')

    def init():
        line.set_data([], [])
        return line

    def make2DPlot(solverData, solver, time, group):
        timeStepData = solverData[time]
        x = timeStepData[:,0]
        y = timeStepData[:,1]
        xUnique = np.unique(x)
        yUnique = np.flip(np.unique(y))
        X, Y = np.meshgrid(xUnique,yUnique)
        g = timeStepData[:,group+2]
        TwoDdata = convertTo2Dgrid(50, 50, g)
        oneDdata = TwoDdata[:,25]
        ax.set_title('Group: '+str(group+1) + ' Time: '+str(time))
        line.set_data(yUnique, oneDdata)
        return line

    def animate(i):
        solver = "hyperbolic"
        solverData = data[solver]
        group = 1
        time = list(solverData.keys())[i]
        return make2DPlot(solverData, solver, time, group+1)

    frames = [i for i in range(600)]
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim = animation.FuncAnimation(fig, animate, init_func = init, 
        frames=frames, interval=100)
    anim.save('test2.mp4', progress_callback=lambda i, 
        n: print(f'Saving frame {i} of {n}'))


if __name__ == "__main__":
    plotSingleCellDepletion("singleCellDepletion.out")
