"""
Author: Zack Taylor

Class that will handle data anlytics primarally for post processing
matlab and libowski solutions to compute the errors.
"""

import numpy as np

class computeErrors:

    @staticmethod
    def computeRMSE(anaSolution, approxSolution):
        """
        Computes the RMSE between the matlab solution and libowski 
        solution. Each column in the input arrays is the solution
        at a single time step. each row in the column is a species 
        concentration. The returned array is a vector where each 
        entry is the RMSE of the difference between the matlab 
        solution and the libowski solution at a single time step.
        """
        assert(anaSolution.shape == approxSolution.shape)
        RMSE = []
        # Loop though columns
        for col in range(anaSolution.shape[1]):
            # Loop though spec concentrations
            sum_ = 0.0
            for specID in range(anaSolution.shape[0]):
                ana = anaSolution[specID, col]
                approx = approxSolution[specID, col]
                if ana != 0.0:
                    relativeError = (ana - approx)/ana
                    if relativeError != 1.0:
                        sum_ += relativeError**2.
            RMSEstep = (sum_/anaSolution.shape[0])**0.5
            RMSE.append(RMSEstep)
        return np.asarray(RMSE)

    @staticmethod
    def computeMaxRelError(anaSolution, approxSolution):
        assert(anaSolution.shape == approxSolution.shape)
        maxRelError = []
        # Loop though columns
        for col in range(anaSolution.shape[1]):
            # Loop though spec concentrations
            max_ = 0.0
            for specID in range(anaSolution.shape[0]):
                ana = anaSolution[specID, col]
                approx = approxSolution[specID, col]
                if ana != 0.0:
                    relativeError = abs(ana-approx)/ana
                    if relativeError != 1.0:
                        max_ = max(relativeError, max_)
            maxRelError.append(max_)
        return np.asarray(maxRelError)
