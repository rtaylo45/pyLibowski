from mpl_toolkits import mplot3d
import matplotlib
matplotlib.use("Agg")
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
#import seaborn as sns; sns.set()
from dataReaders.reader import read2Ddata, convertTo2Dgrid
from analytics.pyAnalytics import computeErrors
from plottingTools.plotter import make2DGIF
import matplotlib.animation as animation
fname = sys.argv[1]
data = read2Ddata(fname)

def getGroupMaxValues(data):
    solver = list(data.keys())[0]
    solverData = data[solver]
    gMaxs = [-1e10, -1e10, -1e10, -1e10, -1e10, -1e10]
    for time in solverData.keys():
        timeStepData = solverData[time]
        for group in range(2, 8):
            g = timeStepData[:,group]
            for groupValue in g:
                #print(time, group-2, gMaxs[group-2], groupValue)
                gMaxs[group-2] = max(gMaxs[group-2], groupValue) 

    return gMaxs

maxVals = getGroupMaxValues(data)

def animate(i):
    solver = "hyperbolic"
    solverData = data[solver]
    group = 1
    time = list(solverData.keys())[i]
    #return makeHeatMap(solverData, solver, time, group+1)
    #return make3DPlot(solverData, solver, time, group+1)
    return make2DPlot(solverData, solver, time, group+1)

def makeHeatMap(solverData, solver, time, group):
    timeStepData = solverData[time]
    x = timeStepData[:,0]
    y = timeStepData[:,1]
    xUnique = np.unique(x)
    yUnique = np.flip(np.unique(y))
    X, Y = np.meshgrid(xUnique,yUnique)
    xticks = np.empty(5, dtype=object)
    yticks = np.empty(20, dtype=object)
    #xticks[0], xticks[24], xticks[49] = '0', '25', '50'
    g = timeStepData[:,group]
    TwoDdata = convertTo2Dgrid(5, 20, g)
    heatmap = sns.heatmap(TwoDdata, cmap='viridis', xticklabels=xticks,
        yticklabels=yticks, cbar=False, vmin=0, vmax=maxVals[group-2])
    plt.title(solver+" time = "+time+" Group: "+str(group-1))
    return heatmap

def make3DPlot(solverData, solver, time, group):
    timeStepData = solverData[time]
    x = timeStepData[:,0]
    y = timeStepData[:,1]
    xUnique = np.unique(x)
    yUnique = np.flip(np.unique(y))
    X, Y = np.meshgrid(xUnique,yUnique)
    g = timeStepData[:,group]
    ax = plt.axes(projection='3d')
    TwoDdata = convertTo2Dgrid(50, 50, g)
    ax.set_xlabel('x')
    ax.set_ylabel('y') 
    ax.set_zlabel('z');
    ax.set_title("Group: " + str(group-2))
    surface = ax.plot_surface(X, Y, TwoDdata, rstride=1, cstride=1,
        cmap='viridis', edgecolor='none', vmin=0, vmax=maxVals[group-2])
    ax.set_zlim(0, maxVals[group-2])
    ax.set_ylim(0,400)
    ax.set_xlim(0,50)
    ax.view_init(azim=-30)
    plt.title(solver+" time = "+time+" Group: "+str(group-1))
    return ax



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
            xticks = np.empty(50, dtype=object)
            yticks = np.empty(100, dtype=object)
            xticks[0], xticks[24], xticks[49] = '0', '25', '50'
            yticks[0], yticks[49], yticks[99] = '400', '200', '0'
            fig, axs = plt.subplots(2,3, constrained_layout=True, figsize = (12,8))
            ax1, ax2, ax3 = axs[0,0], axs[0,1], axs[0,2]
            ax4, ax5, ax6 = axs[1,0], axs[1,1], axs[1,2]
            axsflat = [ax1, ax2, ax3, ax4, ax5, ax6]
            for group in range(6):
                ax = axsflat[group]
                g = timeStepData[:,group+2]
                g = g/np.sum(g)
                TwoDdata = convertTo2Dgrid(50, 100, g)
                #fig = plt.figure()
                #ax = plt.axes(projection='3d')
                #ax.plot_surface(X, Y, TwoDdata, rstride=1, cstride=1,
                #    cmap='viridis', edgecolor='none')
                #ax.set_xlabel('x')
                #ax.set_ylabel('y') 
                #ax.set_zlabel('z');
                ax.set_title("Group: " + str(group+1))
                heatmap = sns.heatmap(TwoDdata, cmap='viridis', xticklabels=xticks,
                    yticklabels=yticks, ax = ax, cbar = False)
                #fig.suptitle(solver+" time = "+time)
                fig.suptitle("Time = 60 (sec)")
                #plt.xlabel('cm')
                #plt.gca().invert_yaxis()
                #plt.savefig('Group'+str(group-1)+'time='+str(time)+'.png')
                #plt.close()
            #plt.show()
            plt.savefig('neutronPrecursors60sec.png')

def plotSingleStepMatlabError():
    cwd = os.getcwd()
    dataFolder = cwd+'/results/caseStudyNeutronPrecorsors/velocity25/10secondruns'
    dataFoldertemp = cwd+'/results/caseStudyNeutronPrecorsors/velocity60/10secondruns'
    fileBaseName = "caseStudyNeutronprecursors"
    substeps = [0, 2, 4, 6, 8, 10, 12]
    solvers = ['CRAM', 'parabolic', 'hyperbolic', 'pade-method1',
        'pade-method2', 'taylor']
    cauchy = ['CRAM', 'parabolic', 'hyperbolic']
    times = [i for i in range(1,61)]
    times = [10, 20, 30, 40, 50, 60]
    dirs = [dataFolder, dataFoldertemp]
    folderNames = ['velocity 25', 'velocity 60']

    # Plots the errors for velocity 25 and 60 for the first time step
    for solver in cauchy:
        for folderIndex, folder in enumerate(dirs):
            os.chdir(folder)
            solutionfname = fileBaseName+"Solution.csv"
            solution = np.genfromtxt(solutionfname, delimiter=',') 
            errors = []
            for substep in substeps:
                fnamecram = fileBaseName+solver+'Substeps'+str(substep)+'.csv'
                approxSol = np.genfromtxt(fnamecram, delimiter=',')
                l1ErrorArray = computeErrors.computeRMSE(solution, approxSol)
                errors.append(l1ErrorArray[0])
            plt.scatter(substeps, errors, label=folderNames[folderIndex])
        plt.grid()
        plt.yscale('log')
        plt.ylabel('RMSE')
        plt.xlabel('Substeps')
        plt.legend()
        plt.title(solver)
        plt.savefig(solver+'velocity25vsvelocity60.png')
        plt.close()
    """
    # change working directroy
    os.chdir(dataFolder)
    solutionfname = fileBaseName+"Solution.csv"
    solution = np.genfromtxt(solutionfname, delimiter=',') 

    # Plots all solvers on one plot for the first time step
    for solver in cauchy:
        errors = []
        for substep in substeps:
            fnamecram = fileBaseName+solver+'Substeps'+str(substep)+'.csv'
            approxSol = np.genfromtxt(fnamecram, delimiter=',')
            l1ErrorArray = computeErrors.computeRMSE(solution, approxSol)
            errors.append(l1ErrorArray[0])
        plt.scatter(substeps, errors, label=solver)
    plt.grid()
    plt.yscale('log')
    plt.ylabel('RMSE')
    plt.xlabel('Substeps')
    plt.legend()
    plt.savefig('velocity60steps10step1.png')
    plt.close()


    # Plots each solver and the reduction in error with substeps
    substeps = [0, 2, 4, 6, 12]
    for solver in cauchy:
        for substep in substeps:
            fnamecram = fileBaseName+solver+'Substeps'+str(substep)+'.csv'
            approxSol = np.genfromtxt(fnamecram, delimiter=',')
            l1ErrorArray = computeErrors.computeRMSE(solution, approxSol)
            plt.scatter(times, l1ErrorArray, label='Substeps = '+str(substep))
        plt.grid()
        plt.yscale('log')
        plt.ylabel('RMSE')
        plt.xlabel('Time (sec)')
        plt.legend()
        plt.title(solver)
        plt.savefig(solver+'velocity60steps10.png')
        plt.close()

    # Plots all solvers on one figure 
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
    plt.ylabel('RMSE')
    plt.xlabel('Time (sec)')
    plt.legend()
    plt.savefig('velocity60steps10.png')
    plt.close()
    """

def makeGIF():
    fig = plt.figure()
    frames = [i for i in range(600)]
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    anim = animation.FuncAnimation(fig, animate, frames=frames,
                              interval=100)
    anim.save('test2.mp4', progress_callback=lambda i, 
        n: print(f'Saving frame {i} of {n}'))

if __name__ == "__main__":
    #groupMaxs = getGroupMaxValues(data)
    #plotHeatMap(data)
    #plotSingleStepMatlabError()
    make2DGIF(data, maxVals)
    #plotSinglePoint(data)

