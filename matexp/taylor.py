import numpy as np
from numpy.linalg import norm
from pandas import read_excel

def normAm(A, m):
    C = A
    for i in range(m):
        C = C @ A
    return norm(C,1)


def paramaters(A, b, **kwargs):
    m_max = kwargs.get('m_max', 55)
    p_max = kwargs.get('p_max', 8)
    shift = kwargs.get('shift', True)
    forceEstm = kwargs.get('forceEstm', False)
    theta = read_excel('theta.xlsx').to_numpy()
    n = A.shape[0]

    if shift:
        mu = A.trace()/n
        A = A - mu*np.identity(n)

    if not forceEstm:
        normA = norm(A,1)

    if not forceEstm and normA <= 4.*theta[m_max]*p_max*(p_max+3)/(m_max*b.shape[0]):
        c = normA
        alpha = c*np.ones((p_max-1,1))
    else:
        mv = 0
        eta = np.zeros((p_max,1))
        alpha = np.zeros((p_max-1,1))
        for p in range(1,p_max+1):
            c = normAm(A, p+1)
            c = c**(1./(p+1))
            eta[p-1] = c

        for p in range(1,p_max):
            alpha[p-1] = max(eta[p-1],eta[p])

    M = np.zeros((m_max, p_max-1))
    for p in range(2, p_max+1):
        for m in range(p*(p-1)-1, m_max+1):
            M[m-1, p-2] = alpha[p-2]/theta[m-1]
    return M, alpha

def taylor(A, t, b, **kwargs):
    shift = kwargs.get('shift', True)
    M = kwargs.get('M', -1)
    fullTerm = kwargs.get('fullTerm', False)
    tol = 2.**(-53)
    n = A.shape[0]

    if shift:
        mu = A.trace()/n
        A = A - mu*np.identity(n)

    if M == -1:
        tt = 1.
        M, alpha = paramaters(A*t, b) 
    else:
        tt = t

    s = 1
    if t == 0:
        m = 0
    else:
        m_max = M.shape[0]
        p = M.shape[1]
        U = np.diag([x for x in range(1,m_max+1)])
        C = (np.ceil(np.abs(tt)*M)).T @ U
        C[C == 0] = np.inf
        # it is assumed that there is only one min value
        cost = np.amin(C)
        # the (row, col) where the cost is
        m = np.where(C == cost)
        # want the column where the cost is [1] gets the column
        # indexes. They are all the same value so just grab the
        # first index [0]
        m = m[1][0] + 1
        if cost == np.inf:
            cost = 0
        s = int(max(cost/m, 1.))

    if shift:
        eta = np.exp(t*mu/s)
    else:
        eta = 1.

    f = b
    for i in range(1, int(s)+1):
        c1 = norm(b, np.inf)
        for k in range(1, m+1):

            b = (t/(s*k))*(A @ b)
            f = f + b
            c2 = norm(b, np.inf)
            print(k, t/(s*k), b)
            if not fullTerm:
                if c1 + c2 <= tol*norm(f, np.inf):
                    break
                c1 = c2
        f = eta*f
        b = f
    return f
       
def expmv(A, t, n0):
    n = taylor(A, t, n0)
    return n
