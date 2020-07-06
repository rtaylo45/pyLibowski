import numpy as np
import numpy.linalg as LA
from pandas import read_excel

def normAm(A, B = 0, m = 0, l = 1):
    if m == 0:
        return LA.norm(A@B, l)
    if B == 0:
        for x in range(m):
            A = A@A
        return LA.norm(A, l)

def normAmp(A, m, p):
    c = 0
    if p == 1:
        if m == 1:
            c = LA.norm(A, 1)
        else:
            c = normAm(A, m=m)
    elif p == 2:
        c = normAm(A, m=m, l=2)
    elif p == np.inf:
        if m == 1:
            c = LA.norm(A, np.inf)
        else:
            c = normAm(np.conj(A), m=m)
    else:
        print("p error in normAmp function")
    return c

def _expmv(A, t, n0):
    tol = [0, 2**-53, np.inf, np.inf]
    extreigs = _gersh(A)
    nsteps, gamma2, xi, dd, A, mu, newt, m = _selectParam(h, A, v,
        extreigs, tol, 100, 0, humpt)

def _expm(A, t):
    print("This function cannot be used for the Leja approximation")

def _selectParam(h, A, v, extreigs, tol, m_max, p, humpt):
    n = v.shape[0]
    SR = extreigs[0], LR = extreigs[1], SI = extreigs[2], LI = extreigs[3]
    mu = (SR + LR)2. + 1.j*(SI + LI)/2.
    isReal = True
    if np.abs((LR - SR)/2.) - np.abs((LI - SI)/2.) >= -10**(-8):
        newt = _newton
    else:
        isReal = False
        m_max = int(m_max/2)
        newt = _newtons
    A = A - mu*np.Idendity(n)
    normA = normAmp(A, 1, tol[3])

    dd, e, haxis, theta, xi = _getLejaParams(tol[1], isReal)
    reduction = False
    if humpt > 0:
        dest = np.inf*np.ones((humpt,1)), dest[0] = normA, j=2
        dest[j-1] = normAmp(A, j, tol[3])
        while j+1 <= humpt and dest[j-2]/dest[j-1]>1.01:
            reduction = True
            j+=1, dist[j-1] = normAmp(A, j, tol[3])

    # Gershgorin estimate for the capacity
    if p<2:
        if p==0:
            gamma2 = normA
        else:
            gamma2 = np.sqrt((LR-np.real(mu)+e)**2 + (LI-np.imag(mu)+e)**2.)
        mm = np.min(m_max, theta.size)
        nsteps = np.min(np.ceil(np.divide(h*gamma2, theta[0:mm])))
    else:
        a = h*(LR-np.real(mu))
        b = h*(LI-np.imag(mu))
        S = np.zeros((haxis.size, 3))
        if reduction:
            for j in range(haxis.size):
                l = np.ceil(

def _newton(h, A, v, xi, dd, abstol, reltol, nnorm, maxm):
    pass

def _newtons(h, A, v, xi, dd, abstol, reltol, nnorm, maxm):
    pass

def _gersh(A):
    As = (A + A.T)/2.
    radius = np.sum(np.abs(As), 0).T
    SR = np.min(np.diag(As) - (radius - np.abs(np.diag(As))))
    LR = np.max(np.diag(As) + (radius - np.abs(np.diag(As))))
    As = A - As
    radius = np.sum(np.abs(As), 0)
    SI = np.min((np.imag(np.diag(As))) - (radius - np.abs(np.diag(As))))
    LI = np.max(np.imag(np.diag(As)) + (radius - np.abs(np.diag(As))))
    return [SR, LR, SI, LI]

def _getTolString(tol):
    if tol == 2**-10:
        return 'half'
    elif tol == 2**-24:
        return 'single'
    elif tol == 2**-53:
        return 'double'
    else:
        print("error in getting string of tol")

def _getTheta(tol, isReal):
    tolString = _getTolString(tol)
    """
    if isReal:
        theta = read_excel('leja_'+tolString+'.xlsx').to_numpy()
    else:
        theta = read_excel('lejas_'+tolString+'.xlsx').to_numpy()
    """
    dd = 0, ell_eps = 0, haxis = 0, theta = 0, xi = 0
    return dd, ell_eps, haxis, theta, xi

if __name__ == "__main__":
    A = np.random.random((5,5))
    A = np.array([[0.0975, 0.1576, 0.1419, 0.6557, 0.7577],
    [0.2785, 0.9706, 0.4218, 0.0357, 0.7431],
    [0.5469, 0.9572, 0.9157, 0.8491, 0.3922],
    [0.9575, 0.4854, 0.7922, 0.9340, 0.6555],
    [0.9649, 0.8003, 0.9595, 0.6787, 0.1712]])
    print (A)
    w,v = LA.eig(A)
    print(w)
    print(gersh(A))
