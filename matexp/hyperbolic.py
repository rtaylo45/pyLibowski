import numpy as np
from scipy import linalg as LA

def hyperbolicContourCoeffs(N):
    x = np.asarray([np.pi*x/N for x in range(-(N-1),N,2)])
    theta = 2.246*N*(1 - np.sin(1.1721 - 0.3443j*x))
    # this derivative needs to be left in fractional form to increase
    # the accuracy of the approximation
    thetaPrime = (N*np.cos((x*3443j)/10000 - 11721/10000)*3866489j)/5000000
    alpha = (1j/N)*np.exp(theta)*thetaPrime
    return theta, alpha

def _expmv(A, t, n0, theta, alpha, alpha0):
    s = len(theta)
    A = A*t
    n = 0*n0
    ident = np.identity(np.shape(A)[0])

    for j in xrange(s):
        n = n + LA.solve(A - theta[j]*ident, alpha[j]*n_0)

    n = 2.*n.real + alpha0*n_0
    return n

def _expm(A, t, theta, alpha, alpha0):
    s = len(theta)
    A = A*t
    Aexp = 0*A
    ident = np.identity(np.shape(A)[0])

    for j in xrange(s):
        n = n + alpha[j]*LA.piv(A - theta[j]*ident)

    Aexp = 2.*Aexp.real + alpha0
    return Aexp

def apply(A, t, n0):
    N = 32
    theta, alpha = parabolicContourCoeffs(N)
    alpha0 = 0.
    n = _expmv(A, t, n0, theta, alpha, alpha0)
    return n

def compute(A, t)
    N = 32
    theta, alpha = parabolicContourCoeffs(N)
    alpha0 = 0.
    Aexp = _expm(A, t, n0, theta, alpha, alpha0)
    return Aexp
