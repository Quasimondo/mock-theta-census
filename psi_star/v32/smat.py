"""The matrix S(a, x) of state integrals in double precision, its reflection S^c, and the three formulas for Omega_gamma (addenda v30 to v32).
   Sany is the function of psi_star/v31/c1_continuation.py."""
import numpy as np
from mpmath import mpc
from common import Setup
import frames
from npint import NP
from lowint import Low


def tonp(M): return np.array([[complex(M[i, j]) for j in range(M.cols)] for i in range(M.rows)])

def Cmat(a, x):
    p = np.exp(-1j*np.pi/x); T = np.exp(1j*np.pi*(2*a - 1)/x)
    return np.array([[1, -1, 0], [0, p, -p/T], [0, -1, 0]])

def Sraw(a, x):
    um = max(abs(a/np.sqrt(x)), abs((a - 1)/np.sqrt(x)), abs((a + x)/np.sqrt(x)), abs((a + x - 1)/np.sqrt(x))) + 0.2
    if x.imag < 0: return Low(x, umax=um).Smat(a)
    assert x.real > 0, "npint.py needs Re x > 0 when Im x >= 0"
    d = min(0.75, 0.9*(np.pi/2 - abs(np.angle(np.sqrt(x)))))
    g = np.pi*np.sin(2*d); R = (2*np.pi*um + np.sqrt((2*np.pi*um)**2 + 4*g*36.0))/(2*g)
    return NP(x, delta=d, R=R).Smat(a)

def Sany(a, x):
    """S(a, x) for x off the negative real axis (Im x >= 0 needs Re x > 0); Re a is first reduced to [0, 1) with the dual shift."""
    x = complex(x); a = complex(a)
    k = int(np.floor(a.real)); a0 = a - k
    S = Sraw(a0, x)
    for j in range(k): S = Cmat(a0 + j, x) @ S
    for j in range(-k): S = np.linalg.solve(Cmat(a0 - j - 1, x), S)
    return S

def Sc(a, x):
    """S^c(a, x) = conj S(-conj a, -conj x)"""
    return np.conj(Sany(-np.conj(a), -np.conj(x)))

def gam(a, x): return a/(1 + 2*x), x/(1 + 2*x)

def lam_np(a, x):
    ap, xp = gam(a, x)
    return np.exp(1j*np.pi/12 - 5j*np.pi*x*xp/6 - 1j*np.pi*(a - ap) - 1j*np.pi*a*ap)/np.sqrt(1 + 2*x)

def Om_series(a, x):
    """Im x > 0: J(a;x)^-1 Lambda^-1 J(a';x') by series (mpmath)"""
    ap, xp = gam(a, x); l = lam_np(a, x)
    J = tonp(frames.frame_q(Setup(mpc(x)), mpc(a))); Jp = tonp(frames.frame_q(Setup(mpc(xp)), mpc(ap)))
    return np.linalg.solve(J, np.diag([1, 1/l, -1/l]) @ Jp)

def Om1(a, x):
    """S(a,x)^-1 S(a',x'):  x and x' off (-inf, 0]"""
    ap, xp = gam(a, x); return np.linalg.solve(Sany(a, x), Sany(ap, xp))

def Om2(a, x):
    """S^c(a,x)^-1 S^c(a',x'):  x and x' off [0, inf)"""
    ap, xp = gam(a, x); return np.linalg.solve(Sc(a, x), Sc(ap, xp))

def Om3(a, x):
    """S(a, x+1)^-1 S^c(a', x'-1):  x off (-inf, -1/2]"""
    ap, xp = gam(a, x); return np.linalg.solve(Sany(a, x + 1), Sc(ap, xp - 1))
