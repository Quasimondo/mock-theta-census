"""Faddeev's quantum dilogarithm near real b, in the variable w = 2 pi b y (so y = b s, w = i t s = 2 pi x s; x = b^2, t = -2 pi i x).
   Lam(w) = log Phi_b(w/(2 pi b)) = int_{R + i0} e^(-i w sigma) / (4 sigma sinh(pi sigma) sinh(pi x sigma)) d sigma,   |Im w| < pi (1 + Re x), Re x > 0.
Outside the strip: functional equation in the direction i/b (w -> w + 2 pi i), with dual variable w' = e^(w/x) and p = e^(-pi i/x):
   Phi_b(y) = Phi_b(y - iN/b) / prod_{j<N} (1 + p^(2j+1) w'),      Phi_b(y) = Phi_b(y + iM/b) prod_{j<M} (1 + p^(-(2j+1)) w').
Right half-plane: inversion  Lam(w) = -t/24 + pi^2/(6t) + w^2/(2t) - Lam(-w).
Two implementations: mpmath (high precision, in the strip) and numpy (double precision, any w)."""
import numpy as np
from mpmath import mp, mpf, mpc, exp, sinh, pi, quad, qp, sqrt, log, polylog
I = mpc(0, 1)

def Lam_mp(w, x, c=mpf(1)/2):
    """mpmath, |Im w| <= pi, on the line Im sigma = c (0 < c < 1)."""
    f = lambda a: exp(-I*w*(a + I*c))/(4*(a + I*c)*sinh(pi*(a + I*c))*sinh(pi*x*(a + I*c)))
    rate = pi*(1 + x.real) - abs(w.imag); A = (mp.dps*mpf("2.31") + 10)/rate
    pts = [-A] + [-A*mpf(k)/8 for k in range(7, 0, -1)] + [-2, -1, 0, 1, 2] + [A*mpf(k)/8 for k in range(1, 8)] + [A]
    return quad(f, sorted(set(pts)))

def Phi_prod_mp(w, x):
    """product formula (Im x != 0), mpmath; argument w = 2 pi b y."""
    t = -2*pi*I*x; wd = exp(w/x)
    if x.imag > 0:
        q = exp(-t); p = exp(-pi*I/x)
        return qp(-exp(-t/2 + w), q, maxterms=10**6)/qp(-p*wd, p*p)
    q1 = exp(t); p1 = exp(pi*I/x)
    return qp(-p1*wd, p1*p1)/qp(-exp(t/2 + w), q1, maxterms=10**6)

_GRID = {}
def _grid(x, c=0.5, h=0.08):
    """nodes sigma = a + i c and the w-independent part of the integrand, written without overflow:
       1/(sinh(pi s) sinh(pi x s)) = 4 e^(-sg pi (1+x) s) / ((1 - e^(-2 sg pi s)) (1 - e^(-2 sg pi x s))),  sg = sign(Re s)."""
    key = (complex(x), c, h)
    if key not in _GRID:
        A = 40.0/(np.pi*x.real) + 20
        a = np.arange(-A, A + h/2, h); sig = a + 1j*c; sg = np.where(a >= 0, 1.0, -1.0)
        dec = -sg*np.pi*(1 + x)*sig
        fac = h/(sig*(1 - np.exp(-2*sg*np.pi*sig))*(1 - np.exp(-2*sg*np.pi*x*sig)))
        _GRID[key] = (sig, dec, fac)
    return _GRID[key]

def Lam_np(w, x):
    """numpy double precision, Re w <= 0, |Im w| <= pi; trapezoid on Im sigma = c with step 0.16 c, c = min(1, Re(1/x))/2 (v30; before: c = 1/2, step 0.08)."""
    # v30: the line Im sigma = c must stay below the first zeros i and i/x of the two sinh factors
    c = 0.5*min(1.0, (1/complex(x)).real)
    sig, dec, fac = _grid(x, c=c, h=0.16*c)
    return np.sum(np.exp(-1j*w*sig + dec)*fac)

def logPhi_np(w, x):
    """log Phi_b for any w (numpy double), Re x > 0."""
    w = complex(w); x = complex(x); t = -2j*np.pi*x
    if w.real > 0:
        return -t/24 + np.pi**2/(6*t) + w*w/(2*t) - logPhi_np(-w, x)
    N = int(np.round(w.imag/(2*np.pi))); w2 = w - 2j*np.pi*N
    out = Lam_np(w2, x)
    if N != 0:
        lw = w/x; lp = -1j*np.pi/x               # log w', log p
        for j in range(abs(N)):
            e = (2*j + 1)*lp*(1 if N > 0 else -1) + lw
            term = np.log1p(np.exp(e)) if e.real < 300 else e
            out += -term if N > 0 else term
    return out
