import numpy as np
from scipy import linalg as LA

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
    Anrom = norm(A,1)
    s = 0
    A2 = A@A
    I = np.identity(A.shape[0])

    if Anorm < 1.495585217958292e-002:
        u,v = _pade3(A, I, A2)
    elif Anorm < 2.539398330063230e-001:
        A4 = A2@A2
        u,v = _pade5(A, I, A2, A4)
    elif Anrom < 9.504178996162932e-001:
        A4 = A2@A2
        A6 = A4@A2
        u,v = _pade7(A, I, A2, A4, A6)
    elif Anorm < 2.097847961257068e+000:
        A4 = A2@A2
        A6 = A4@A2
        A8 = A6@A2
        u,v = _pade9(A, I, A2, A4, A6, A8)
    else:
        maxnorm = 5.371920351148152
        s = max(0, int(ceil(log2(Anorm/ maxnorm))))
        A = A / 2**s
        A2 = A@A
        A4 = A2@A2
        A6 = A4@A2
        u,v = _pade13(A, I, A2, A4, A6)

    p = u + v  # p_m(a) : numerator
    q = -u + v # q_m(a) : denominator
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
