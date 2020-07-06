import numpy as np
from scipy import linalg as LA
from scipy.special import binom
from scipy.misc import factorial

def normAm(A, B = 0, m = 0):
    if m == 0:
        return LA.norm(A@B, 1)
    if B == 0:
        for x in range(m):
            A = A@A
        return LA.norm(A, 1)

def ell(A, m):
    p = 2*m + 1

    c = np.abs(1./(binom(2*p, p)*factorial(2*p+1)))
    u = 2**-53
    A1normMatrixPower = normAm(np.abs(A), p)
    A1norm = l1norm(A)
    alpha = c*A1normMatrixPower/A1norm
    log2AlphaOverU = np.log2(alpha/u)
    value = int(np.ceil(log2AlphaOverU/(2*m)))
    return np.max(value, 0)

def _pade3(A, I, A2):
    b = [120., 60., 12., 1.]
    temp = b[3]*A2 + b[1]*I
    u = A @ temp
    v = b[2]*A2 + b[0]*I
    return u,v

def _pade5(A, I, A2, A4):
    b = [30240., 15120., 3360., 420., 30., 1.]
    temp = b[5]*A4 + b[3]*A2 + b[1]*I
    u = A @ temp
    v = b[4]*A4 + b[2]*A2 + b[0]*I
    return u,v

def _pade7(A, I, A2, A4, A6):
    b = [17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.]
    temp = b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*I
    u = A @ temp
    v = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*I
    return u,v

def _pade9(A, I, A2, A4, A6, A8):
    b = [17643225600., 8821612800., 2075673600., 302702400., 30270240.,
                2162160., 110880., 3960., 90., 1.]
    temp = b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*I
    u = A @ temp
    v = b[8]*A8 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*I
    return u,v

def _pade13(A, I , A2, A4, A6):
    b = [64764752532480000., 32382376266240000., 7771770303897600.,
    1187353796428800., 129060195264000., 10559470521600., 670442572800.,
    33522128640., 1323241920., 40840800., 960960., 16380., 182., 1.]
    v = b[13]*A6 + b[11]*A4 + b[9]*A2
    temp = A6 @ v
    temp += b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*I
    u = A @ temp
    v = A6 @ temp
    v = += b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*I
    return u,v

def c(order, k):
    k = float(k)
    order = float(order)
    if k == 0:
        return 1
    else:
        return c(order,k-1)*(order+1.-k)/((2.*order+1.-k)*k)

def _expmv(A, t, n0):
    n = _expm(A, t) @ n0
    return n

def _expm(A, t):
    A = A*t
    alpha = 0
    I = np.Identity(A.shape[0])
    theta13 = 4.25

    # try pade3
    A2 = A@A
    d6 = normAm(A2, m=3)**1./6.
    eta1 = np.max(normAm(A2, m=2), d6)
    if (eta1 < 1.495585217958292e-002 and ell(A, 3) == 0):
        U, V = _pade3(A, I, A2)
        p = U + V  # p_m(a) : numerator
        q = -U + V # q_m(a) : denominator
        R = LA.solve(q,p)
        return R

    #try pade5
    A4 = A2@A2
    d4 = np.norm(A4, 1)**1./4.
    eta2 = np.max(d4, d6)
    if (eta < 2.539398330063230e-001 and ell(A, 5) == 0):
        U, V = _pade5(A, I, A2, A4)
        p = U + V  # p_m(a) : numerator
        q = -U + V # q_m(a) : denominator
        R = LA.solve(q,p)
        return R

    # try pade 7 and 9
    A6 = A4@A2
    d6 = l1norm(A6)**1./6.
    d8 = normAm(A4, m=2)**1./8.
    eta = np.max(d6, d8)
    if (eta3 < 9.504178996162932e-001 adn ell(A,7) == 0):
        U, V = _pade7(A, I, A2, A4, A6)
        p = U + V  # p_m(a) : numerator
        q = -U + V # q_m(a) : denominator
        R = LA.solve(q,p)
        return R

    if (eta3 < 2.097847961257068e+000 and ell(A,9) == 0):
        A8 = A6@A2
        U, V = _pade9(A, I, A2, A4, A6, A8)
        p = U + V  # p_m(a) : numerator
        q = -U + V # q_m(a) : denominator
        R = LA.solve(q,p)
        return R

    # Do pade 13
    d10 = normAm(A4,B=A6)
    eta4 = np.max(d8,d10)
    eta5 = np.min(eta3, eta4)
    log2EtaOverTheta = np.log2(eta5/theta13)
    value = int(np.ceil(log2EtaOverTheta))
    # find number of squarings
    alpha = np.max(value, 0)
    alpha = alpha + ell(A*(2**-alpha), 13)
    # scale the matrix
    B = A*(2**-alpha)
    B2 = A2*(2**-alpha)
    B4 = A4*(2**-alpha)
    B6 = A6*(2**-alpha)
    U, V = _pade13(B, I, B2, B4, B6)

    p = U + V  # p_m(a) : numerator
    q = -U + V # q_m(a) : denominator
    R = LA.solve(q,p)
    # squaring step to undo scaling
    for i in range(s):
        R = dot(R,R)

    return R

def apply(A, t, n0):
    R = _expmv(A, t, n0)
    return R

def compute(A, t):
    Aexp = _expm(A, t)
    return Aexp
