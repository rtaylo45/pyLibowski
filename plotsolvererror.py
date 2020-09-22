import numpy as np
from analytics.pyAnalytics import computeErrors
import matplotlib.pyplot as plt

# Base file name
basename = "results/casestudy1/caseStudy1"
#basename = "results/mole/moleProblem12"

# Solver names
solvers = ["CRAM", "parabolic", "hyperbolic", "pade-method1", "pade-method2", "taylor"]

solution = np.genfromtxt(basename+"Solution.csv", delimiter=',', dtype=np.float64)
# Plots the RMSE for all the solvers
fig = plt.figure()
for solver in solvers:
    fname = basename + solver + ".csv"
    approx = np.genfromtxt(fname, delimiter=',', dtype=np.float64)
    error = computeErrors.computeRMSE(solution, approx)
    time = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
    plt.scatter(time, error, label=solver, s=30.0)
plt.legend(bbox_to_anchor=(0.5, 1.00), loc="upper center",
                bbox_transform=fig.transFigure, ncol=3)
plt.grid()
plt.ylabel("RMSE")
plt.xlabel("Time (days)")
#lt.title("Mole Problem 12")
#lt.title("Case Study One")
#plt.ylim(10**-17, 10**2)
plt.yscale("log")
#plt.savefig("molePRMSE.png",dpi=500)
#plt.close()
plt.show()

"""
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
#plt.savefig("caseStudy1MaxRelativeError.png",dpi=500)
#plt.close()
plt.show()


# Plots the RMSE for the Cauchy solvers with substeps
substeps = ["0", "2", "4", "6", "12"]
solvers = ["CRAM", "parabolic", "hyperbolic"]
for solver in solvers:
    #fig = plt.figure(figsize=(12,8))
    fig = plt.figure()
    for substep in substeps:
        fname = basename + solver + "Steps" + substep + ".csv"
        approx = np.genfromtxt(fname, delimiter=',', dtype=np.float64)
        error = computeErrors.computel1RelError(solution, approx)
        time = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
        plt.scatter(time, error, label="Substeps = "+substep, s=20.0)
    plt.legend(bbox_to_anchor=(0.5, 1.00), loc="upper center",
                bbox_transform=fig.transFigure, ncol=3)
    plt.grid()
    plt.ylabel(r"$l_2$ Relative Error", fontsize=15)
    plt.xlabel("Time (days)", fontsize=15)
    #plt.title("Mole Problem 12")
    #plt.title(solver)
    #plt.title("Case Study One "+solver+" Solver")
    #plt.ylim(10**-17, 10**2)
    plt.yscale("log")
    #plt.savefig("caseStudy1RMSE"+solver+"Substeps"+substep+".png",dpi=500)
    plt.savefig("decay"+solver+"Substeps"+substep+".png",dpi=500)
    plt.close()
    #plt.show()
"""

# Plots the error for each of the cauchy solvers at the first time stepj
substeps = ["0", "2", "4", "6", "8", "10", "12"]
solvers = ["CRAM", "parabolic", "hyperbolic"]
fig = plt.figure()
for solver in solvers:
    errors = []
    #fig = plt.figure(figsize=(12,8))
    for substep in substeps:
        fname = basename + solver + "Steps" + substep + ".csv"
        approx = np.genfromtxt(fname, delimiter=',', dtype=np.float64)
        error = computeErrors.computeRMSE(solution, approx)
        errors.append(error[0])
    plt.scatter([int(step) for step in substeps], errors, label=solver, s=30.0)
#plt.legend(bbox_to_anchor=(0.5, 0.95), loc="upper center",
                #bbox_transform=fig.transFigure, ncol=3)
plt.legend()
#plt.ylim(10**-17, 10**-6)
plt.grid()
plt.ylabel("RMSE")
plt.xlabel("Substeps")
plt.yscale("log")
#plt.savefig("moleP12CauchySubstepsRMSE.png",dpi=500)
#plt.close()
#plt.show()
