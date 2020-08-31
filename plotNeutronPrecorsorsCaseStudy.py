from mpl_toolkits import mplot3d
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import seaborn as sns; sns.set()
from dataReaders.reader import read2Ddata, convertTo2Dgrid
from analytics.pyAnalytics import computeErrors
import matplotlib.animation as animation

def animate(i):
    fname = "caseStudyNeutronPrecursors.out"
    data = read2Ddata(fname)
    solver = "hyperbolic"
    solverData = data[solver]
    group = 3
    time = list(solverData.keys())[i]
    return makeHeatMap(solverData, solver, time, group+1)

def makeHeatMap(solverData, solver, time, group):
    timeStepData = solverData[time]
    x = timeStepData[:,0]
    y = timeStepData[:,1]
    xUnique = np.unique(x)
    yUnique = np.flip(np.unique(y))
    X, Y = np.meshgrid(xUnique,yUnique)
    xticks = np.empty(25, dtype=object)
    yticks = np.empty(50, dtype=object)
    #xticks[0], xticks[24], xticks[49] = '0', '25', '50'
    g = timeStepData[:,group]
    TwoDdata = convertTo2Dgrid(25, 50, g)
    heatmap = sns.heatmap(TwoDdata, cmap='viridis', xticklabels=xticks,
        yticklabels=yticks)
    plt.title(solver+" time = "+time+" Group: "+str(group-1))
    return heatmap

def plotHeatMap(data):
    os.chdir('../')
    os.chdir(os.getcwd() + '/PhDScripts/results/caseStudyNeutronPrecorsors')
    for solver in data.keys():
        solverData = data[solver]
        for time in solverData.keys():
            timeStepData = solverData[time]
            x = timeStepData[:,0]
            y = timeStepData[:,1]
            xUnique = np.unique(x)
            yUnique = np.flip(np.unique(y))
            X, Y = np.meshgrid(xUnique,yUnique)
            xticks = np.empty(5, dtype=object)
            yticks = np.empty(20, dtype=object)
            #xticks[0], xticks[24], xticks[49] = '0', '25', '50'
            for group in range(4, 5):
                g = timeStepData[:,group]
                TwoDdata = convertTo2Dgrid(5, 20, g)
                #fig = plt.figure()
                #ax = plt.axes(projection='3d')
                #ax.plot_surface(X, Y, TwoDdata, rstride=1, cstride=1,
                #    cmap='viridis', edgecolor='none')
                #ax.set_xlabel('x')
                #ax.set_ylabel('y') 
                #ax.set_zlabel('z');
                heatmap = sns.heatmap(TwoDdata, cmap='viridis', xticklabels=xticks,
                    yticklabels=yticks)
                plt.title(solver+" time = "+time+" Group: "+str(group-1))
                #plt.xlabel('cm')
                #plt.gca().invert_yaxis()
                plt.savefig('Group'+str(group-1)+'time='+str(time)+'.png')
                plt.close()

def plotSingleStepMatlabError(dataFolder):
    fileBaseName = dataFolder+"/caseStudyNeutronprecursors"
    solutionfname = fileBaseName+"Solution.csv"
    solution = np.genfromtxt(solutionfname, delimiter=',') 
    substeps = [0, 2, 4, 6, 12]
    solvers = ['CRAM', 'parabolic', 'hyperbolic', 'pade-method1',
        'pade-method2', 'taylor']
    cauchy = ['CRAM', 'parabolic', 'hyperbolic']

    times = [0.1, 0.5, 1.0, 10.0, 20.0]

    for substep in substeps:
        fnamecram = fileBaseName+'CRAMSubsteps'+str(substep)+'.csv'
        approxSol = np.genfromtxt(fnamecram, delimiter=',')
        l1ErrorArray = computeErrors.computeRMSE(solution, approxSol)
        plt.scatter(times, l1ErrorArray, label='Substeps = '+str(substep))
    plt.grid()
    plt.yscale('log')
    plt.legend()
    plt.title('CRAM')
    plt.show()
    plt.close()

    for solver in solvers:
        substep = '6'
        if solver in cauchy:
            fname = fileBaseName+solver+'Substeps'+substep+'.csv'
        else:
            fname = fileBaseName+solver+'.csv'
        approxSol = np.genfromtxt(fname, delimiter=',')
        l1ErrorArray = computeErrors.computeRMSE(solution, approxSol)
        plt.scatter(times, l1ErrorArray, label=solver)
    plt.grid()
    plt.yscale('log')
    plt.legend()
    plt.title('Single Step Solver Error')
    plt.show()
        

if __name__ == "__main__":
    fig = plt.figure()
    fname = sys.argv[1]
    data = read2Ddata(fname)
    plotHeatMap(data)

    # Set up formatting for the movie files
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    #anim = animation.FuncAnimation(fig, animate, frames=5,
    #                          interval=100, blit = True)
    #anim.save('test.mp4', writer = writer)
