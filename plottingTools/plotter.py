"""
Author: Zack Taylor

Class that lets me plot things in libowski
"""
import matplotlib.pyplot as plt
import numpy as np

class pyPlotter:

    @staticmethod
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
        
    @staticmethod
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



if __name__ == "__main__":
    pyPlotter.plotSingleCellDepletion("singleCellDepletion.out")
