"""Independent implementation of the objects of docs/psi-star-theorems-v1.md
used by the stress tests (double precision, numpy).

Conventions (as in §0.1 of the consolidated document):
  t = |t| e^{i theta},  x = i t / (2 pi),  b = sqrt(x),  q = e^{-t},  p = e^{-2 pi^2 / t}
  w = 2 pi b y = i t s,  w' = e^{2 pi s} = e^{w/x}
  Lambda(w) = log Phi_b(w / (2 pi b))

Two independent evaluations of Lambda:
  lam_prod : product formulas (Im x > 0 and Im x < 0, the latter = Phi_{1/b})
  lam_fad  : Faddeev's integral on Im sigma = 1/2 (trapezoid rule, exponentially
             convergent), reduced to Re w <= 0 by the inversion relation and to
             |Im w| <= pi by the functional equation in the direction i/b.
Both return *a* logarithm (the real part is what the bounds need; exp() is exact).
Nothing here is imported from psi_star/v*.
"""
import numpy as np
from math import pi, log

CL2 = 1.0149416064096536  # Cl_2(pi/3)


def tx(t):
    t = complex(t)
    x = 1j * t / (2 * pi)
    return t, x


def log1pexp(z):
    """log(1 + e^z), principal log of the factor, stable for large Re z."""
    z = np.asarray(z, dtype=complex)
    out = np.empty_like(z)
    big = z.real > 30
    out[big] = z[big] + np.log1p(np.exp(-z[big]))
    out[~big] = np.log1p(np.exp(z[~big]))
    return out


# ---------------------------------------------------------------- products
def lam_prod(w, t, nmax_cap=400000):
    """log Phi_b via the product formulas; valid for Re t != 0.
    Re t > 0 (Im x > 0): Phi = (-q^{1/2} e^w; q)_inf / (-p w'; p^2)_inf
    Re t < 0 (Im x < 0): Phi = (-p1 w'; p1^2)_inf / (-q1^{1/2} e^w; q1)_inf,
                          q1 = e^{t}, p1 = e^{2 pi^2/t}.
    Returns None if more than nmax_cap factors would be needed."""
    t, x = tx(t)
    w = np.atleast_1d(np.asarray(w, dtype=complex))
    sgn = 1.0 if t.real > 0 else -1.0
    t1 = sgn * t                       # Re t1 > 0
    lw2 = w / x                        # log w'
    lp1 = -2 * pi ** 2 / t1            # log p1  (p1 = p for Re t > 0)
    # number of q-factors:  |e^{w-(n+1/2) t1}| < 1e-18
    nq = int(np.ceil((max(w.real.max(), 0) + 42) / t1.real)) + 2
    npf = int(np.ceil((max(lw2.real.max(), 0) + 42) / (-2 * lp1.real) + 1)) + 2
    if nq > nmax_cap or npf > nmax_cap:
        return None
    out = np.zeros(w.shape, dtype=complex)
    n = np.arange(nq)
    chunk = max(1, int(2e7 // max(nq, 1)))
    for i0 in range(0, len(w), chunk):
        ww = w[i0:i0 + chunk, None]
        sq = log1pexp(ww - (n[None, :] + 0.5) * t1).sum(axis=1)
        m = np.arange(npf)
        sp = log1pexp(lw2[i0:i0 + chunk, None] + (2 * m[None, :] + 1) * lp1).sum(axis=1)
        out[i0:i0 + chunk] = sgn * (sq - sp)
    return out


# ---------------------------------------------------------------- Faddeev
def _logsinh(z):
    z = np.asarray(z, dtype=complex)
    out = np.empty_like(z)
    m = np.abs(z.real) < 20
    out[m] = np.log(np.sinh(z[m]))
    p = z.real >= 20
    out[p] = z[p] - log(2) + np.log1p(-np.exp(-2 * z[p]))
    n = z.real <= -20
    out[n] = -z[n] - log(2) + 1j * pi + np.log1p(-np.exp(2 * z[n]))
    return out


class Faddeev:
    """Lambda(w) = int_{R + i/2} e^{-i w sigma} / (4 sigma sinh(pi sigma) sinh(pi x sigma)) d sigma
    for Re w <= 0, |Im w| <= pi; trapezoid rule with step h, truncated per w."""

    def __init__(self, t, h=0.04, c=0.5, tail=46.0):
        self.t, self.x = tx(t)
        self.h, self.c, self.tail = h, c, tail
        x = self.x
        assert x.real > 0
        # poles of 1/sinh(pi x sigma) at i n/x: must stay above the line
        assert (1 / x).real > c + 0.3, "x too large for the line Im sigma = 1/2"
        self.Amax = None
        self._grid(200.0)

    def _grid(self, A):
        self.Amax = A
        a = np.arange(-A, A + self.h / 2, self.h)
        sg = a + 1j * self.c
        x = self.x
        self.a = a
        self.sg = sg
        self.logD = -np.log(4 * sg) - _logsinh(pi * sg) - _logsinh(pi * x * sg)

    def strip(self, w):
        """Lambda for Re w <= 0, |Im w| <= pi (array)."""
        w = np.atleast_1d(np.asarray(w, dtype=complex))
        x = self.x
        out = np.empty(w.shape, dtype=complex)
        gap = pi * (1 + x.real) - np.abs(w.imag)
        # tail of |integrand| ~ e^{-gap |a|} / (4 pi^2 |x| a^2) (for |a| < 1/|x|)
        A = np.minimum((self.tail + max(0.0, log(1.0 / (4 * pi ** 2 * abs(x))))) / gap, 2.0e4)
        if A.max() > self.Amax:
            self._grid(A.max() * 1.05)
        n0 = len(self.a) // 2                      # index of a = 0 (grid symmetric)
        order = np.argsort(A)
        k = 0
        while k < len(order):
            # group points with similar truncation (within 25%)
            Ak = A[order[k]]
            j = k
            while j < len(order) and A[order[j]] <= 1.25 * Ak and j - k < 256:
                j += 1
            idx = order[k:j]
            m = int(A[idx].max() / self.h) + 1
            sl = slice(max(n0 - m, 0), n0 + m + 1)
            E = np.exp(-1j * w[idx, None] * self.sg[None, sl] + self.logD[None, sl])
            out[idx] = self.h * E.sum(axis=1)
            k = j
        return out

    def lam(self, w):
        """Lambda(w) for any w (Re x > 0)."""
        t, x = self.t, self.x
        w = np.atleast_1d(np.asarray(w, dtype=complex))
        out = np.empty(w.shape, dtype=complex)
        neg = w.real > 0
        wl = np.where(neg, -w, w)               # Re wl <= 0
        N = np.round(wl.imag / (2 * pi)).astype(int)
        w0 = wl - 2j * pi * N                  # |Im w0| <= pi
        base = self.strip(w0)
        lp = -2 * pi ** 2 / t                  # log p
        corr = np.zeros(w.shape, dtype=complex)
        for i in np.nonzero(N)[0]:
            lw2 = wl[i] / x                    # log w' of the unshifted point
            n = N[i]
            if n > 0:   # Phi(y) = Phi(y - iN/b) / prod_{j<N} (1 + p^{2j+1} w')
                j = np.arange(n)
                corr[i] = -log1pexp(lw2 + (2 * j + 1) * lp).sum()
            else:       # Phi(y) = Phi(y + iM/b) * prod_{j<M} (1 + p^{-(2j+1)} w')
                j = np.arange(-n)
                corr[i] = log1pexp(lw2 - (2 * j + 1) * lp).sum()
        laml = base + corr
        inv = -t / 24 + pi ** 2 / (6 * t)
        out = np.where(neg, -laml + inv + w ** 2 / (2 * t), laml)
        return out


def lam_any(w, t, fad=None, cap=400000):
    """Lambda by products where cheap, else Faddeev (Re w > 0 via inversion)."""
    t, x = tx(t)
    w = np.atleast_1d(np.asarray(w, dtype=complex))
    neg = w.real > 0
    wl = np.where(neg, -w, w)
    if x.real <= 0:
        cap = 10 ** 8          # Faddeev's w-form needs Re x > 0; products always work there
    lp = lam_prod(wl, t, nmax_cap=cap) if abs(t.real) > 1e-12 else None
    if lp is None:
        if fad is None:
            fad = Faddeev(t)
        lp = fad.lam(wl)
    inv = -t / 24 + pi ** 2 / (6 * t)
    return np.where(neg, -lp + inv + w ** 2 / (2 * t), lp)


# ---------------------------------------------------------------- dilog
from scipy.special import spence as _spence


def li2(z):
    """Principal Li_2 (cut [1, inf)): Li2(z) = spence(1 - z) for complex z."""
    z = np.asarray(z, dtype=complex)
    return _spence(1 - z)


def li2_neg_exp(w):
    """Li_2(-e^w) principal (as a function of the point -e^w)."""
    return li2(-np.exp(np.asarray(w, dtype=complex)))


# ---------------------------------------------------------------- contour, G, psi, H
DEG = pi / 180


def angles(theta_deg, delta_deg=10.0):
    beta = (theta_deg + 90.0) / 2
    return beta, beta - delta_deg, beta + delta_deg   # beta, alpha, gamma (degrees)


def G_of_w(w):
    """G of Part 3 (piecewise)."""
    w = np.atleast_1d(np.asarray(w, dtype=complex))
    out = np.empty_like(w)
    L = w.real <= 0
    out[L] = -w[L] ** 2 / 2 + li2_neg_exp(w[L]) + pi ** 2 / 12
    R = ~L
    out[R] = -w[R] ** 2 - li2_neg_exp(-w[R]) - pi ** 2 / 12
    return out


def psi_of(r, theta_deg, delta_deg=10.0):
    """psi(rho) on the contour, as a function of the signed parameter r (w = |t| r e^{i alpha});
    returned as psi(rho) (independent of |t|).  rho = |r| (|t| = 1 scaling)."""
    beta, al, ga = angles(theta_deg, delta_deg)
    w = np.asarray(r, dtype=float) * np.exp(1j * al * DEG)
    G = G_of_w(w)
    return (np.exp(-1j * theta_deg * DEG) * G).real - pi * np.abs(r) * np.cos(ga * DEG)


def contour_w(rho_signed, theta_deg, delta_deg=10.0):
    beta, al, ga = angles(theta_deg, delta_deg)
    return np.asarray(rho_signed, dtype=float) * np.exp(1j * al * DEG)


def log_cosh_term(s):
    """log(2 cosh(pi s)), stable."""
    s = np.asarray(s, dtype=complex)
    z = pi * s
    sg = np.where(z.real >= 0, 1.0, -1.0)
    zz = sg * z
    return zz + np.log1p(np.exp(-2 * zz))


def logH(rho_signed, tabs, theta_deg, version=1, fad=None, delta_deg=10.0, cap=400000):
    """log H on the contour of Parts 3-5 at the points w = rho_signed e^{i alpha}.
    H = C0 e^{-w^2/(2t)} / (Phi_b(b s) 2 cosh(pi s)),  C0 = 2^{-1/2} e^{pi^2/(12t) + t/24};
    version 0 multiplies by e^{w} (psi_0^* integrand)."""
    t = tabs * np.exp(1j * theta_deg * DEG)
    w = contour_w(rho_signed, theta_deg, delta_deg)
    s = w / (1j * t)
    lam = lam_any(w, t, fad=fad, cap=cap)
    out = -0.5 * log(2) + pi ** 2 / (12 * t) + t / 24 - w ** 2 / (2 * t) - lam - log_cosh_term(s)
    if version == 0:
        out = out + w
    return out, s, w
