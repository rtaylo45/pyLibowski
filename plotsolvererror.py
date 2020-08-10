import numpy as np
from analytics.pyAnalytics import computeErrors
import matplotlib.pyplot as plt

# Base file name
basename = "results/caseStudy1"
# Solver names
solvers = ["CRAM", "parabolic", "hyperbolic", "pade-method1", "pade-method2"]
#solvers = ["CRAM", "parabolic", "hyperbolic", "pade-method1", "pade-method2", "taylor"]
#solvers = ["CRAM"]

solution = np.genfromtxt(basename+"Solution.csv", delimiter=',', dtype=np.float64)

# Plots the RMSE for all the solvers
fig = plt.figure(figsize=(12,8))
for solver in solvers:
    fname = basename + solver + ".csv"
    approx = np.genfromtxt(fname, delimiter=',', dtype=np.float64)
    RMSE = computeErrors.computeRMSE(solution, approx)
    time = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    plt.scatter(time, RMSE, label=solver, s=50.0)
plt.legend(bbox_to_anchor=(1,0), loc="lower right",
                bbox_transform=fig.transFigure, ncol=3)
plt.grid()
plt.ylabel("RMSE")
plt.xlabel("Time (yr)")
#plt.title("Mole Problem 12")
plt.title("Case Study One")
#plt.ylim(10**-17, 10**2)
plt.yscale("log")
plt.savefig("caseStudy1RMSE.png",dpi=500)
plt.close()
#plt.show()

# Plots the max relative error for all the solvers
fig = plt.figure(figsize=(12,8))
for solver in solvers:
    fname = basename + solver + ".csv"
    approx = np.genfromtxt(fname, delimiter=',', dtype=np.float64)
    maxRelError = computeErrors.computeMaxRelError(solution, approx)
    time = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    plt.scatter(time, maxRelError, label=solver, s=50.0)
plt.legend(bbox_to_anchor=(1,0), loc="lower right",
                bbox_transform=fig.transFigure, ncol=3)
plt.grid()
plt.ylabel("Max Relative Error")
plt.xlabel("Time (yr)")
plt.title("Case Study One")
#plt.ylim(10**-17, 10**2)
plt.yscale("log")
plt.savefig("caseStudy1MaxRelativeError.png",dpi=500)
plt.close()
#plt.show()


# Plots the RMSE for the Cauchy solvers with substeps
substeps = ["0", "2", "4", "6", "12"]
solvers = ["CRAM", "parabolic", "hyperbolic"]
for solver in solvers:
    fig = plt.figure(figsize=(12,8))
    for substep in substeps:
        fname = basename + solver + "Substeps" + substep + ".csv"
        approx = np.genfromtxt(fname, delimiter=',', dtype=np.float64)
        RMSE = computeErrors.computeRMSE(solution, approx)
        time = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        plt.scatter(time, RMSE, label="Substeps = "+substep, s=50.0)
    plt.legend(bbox_to_anchor=(1,0), loc="lower right",
                    bbox_transform=fig.transFigure, ncol=3)
    plt.grid()
    plt.ylabel("RMSE")
    plt.xlabel("Time (yr)")
    #plt.title("Mole Problem 12")
    plt.title("Case Study One "+solver+" Solver")
    #plt.ylim(10**-17, 10**2)
    plt.yscale("log")
    plt.savefig("caseStudy1RMSE"+solver+"Substeps"+substep+".png",dpi=500)
    plt.close()
    #plt.show()
