"""Addendum v17, the shift identity and the size of the shifted integral.
With H the integrand of P J (v16 section 4), w = i t s, tau_n(t) = q^(n^2+n)/(2 (-q;q)_n), q = e^(-t):
   P J = sum_{n<N} (-1)^n tau_n(t) + e^(t/8) (-1)^N int_Gamma H(s) D_N(w) ds,
   D_N(w) = q^(N^2) e^(-2Nw)/(-q^(1/2) e^(-w); q)_N = prod_{j=1}^N g(zeta_j),   g(z) = z^2/(1+z),  zeta_j = exp(-(w + t(j - 1/2))).
(psi0 version: tau_n = q^(n^2)/(2 (-q;q)_n), extra factor q^(-N) in D_N.)
 (1) the identity, at real x and on both sides, for several N;
 (2) |P J - sum_{n<N} (-1)^n tau_n| as a function of N at real x: it falls like 2^(-N) and bottoms out near e^(-Cl2(pi/3)/|t|);
 (3) the elementary constants used in the proof.
Double precision (numpy), Phi_b from Faddeev's integral (fad.py).   Usage: python3 g1_shift.py"""
import numpy as np
from fad import logPhi_np
DELTA = 10.0

def logH(r, t, x, case="psi1"):
    beta = np.angle(np.sqrt(x)); ga = beta + np.radians(DELTA)
    s = r*np.exp(-1j*ga); w = 1j*t*s
    z = np.pi*s*(1 if s.real >= 0 else -1)
    out = np.pi**2/(12*t) + t/24 - 0.5*np.log(2) + t*s*s/2 - logPhi_np(w, x) - (z + np.log1p(np.exp(-2*z)))
    return (out + w if case == "psi0" else out), w, ga

def logD(N, w, t, case="psi1"):
    tot = 0j
    for j in range(1, N + 1):
        a = -(w + t*(j - 0.5))                      # log zeta_j
        l1p = (a + np.log1p(np.exp(-a))) if a.real > 0 else np.log1p(np.exp(a))
        tot += 2*a - l1p
    return tot + (N*t if case == "psi0" else 0)

def integral(t, x, N=0, case="psi1", R=22.0, panels=110, order=12):
    xs, ws = np.polynomial.legendre.leggauss(order); tot = 0j
    edges = np.linspace(-R, R, panels + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        for xi, wi in zip(xs, ws):
            r = (a + b)/2 + (b - a)/2*xi
            lh, w, ga = logH(r, t, x, case)
            tot += wi*(b - a)/2*np.exp(lh + (logD(N, w, t, case) if N else 0))*np.exp(-1j*ga)
    return tot*np.exp(t/8 - (t/2 if case == "psi0" else 0))*(-1)**N

def tau(n, t, case="psi1"):
    q = np.exp(-t); den = 2.0
    for j in range(1, n + 1): den = den*(1 + q**j)
    return q**(n*n + (n if case == "psi1" else 0))/den

def tx(T, th):
    t = T*np.exp(1j*np.radians(th)); x = 1j*t/(2*np.pi)
    return t, (complex(x.real, 0.0) if th == -90 else x)

if __name__ == "__main__":
    print("(1) shift identity: |P J - [ sum_{n<N} (-1)^n tau_n + shifted integral ]|")
    for case, T, th in (("psi1", 0.3, -90), ("psi1", 0.1, -90), ("psi1", 0.2, -60), ("psi1", 0.2, -150), ("psi0", 0.2, -90), ("psi0", 0.2, -120)):
        t, x = tx(T, th); pj = integral(t, x, 0, case); out = []
        for N in (1, 2, 5):
            rhs = sum((-1)**n*tau(n, t, case) for n in range(N)) + integral(t, x, N, case)
            out.append(f"N={N}: {abs(pj - rhs):.1e}")
        print(f"    {case}, |t| = {T}, theta = {th:5d}:  P J = {pj:.10f}   " + "   ".join(out))
    print("(2) real x: |P J - sum_{n<N} (-1)^n tau_n(t)| for increasing N")
    from scipy.special import spence
    Cl2 = 1.0149416064096536
    for T in (0.2, 0.1, 0.05):
        t, x = tx(T, -90); pj = integral(t, x, 0); part = 0j; row = []
        Ns = sorted(set([0, 1, 2, 4, 8, 12, 16, 24, 32, 40, 48, int(2*np.pi/(3*T))]))
        n = 0
        for N in range(0, max(Ns) + 1):
            if N in Ns: row.append(f"N={N}: {abs(pj - part):.1e}")
            part += (-1)**N*tau(N, t)
        print(f"    |t| = {T}:  " + "  ".join(row) + f"     e^(-Cl2/|t|) = {np.exp(-Cl2/T):.1e},  2 pi/(3|t|) = {2*np.pi/(3*T):.1f}")
    print("(3) constants")
    ph = np.linspace(0, 2*np.pi, 20001); eta = np.exp(1j*ph)
    print(f"    min over |eta| = 1 of |1 - e^eta| = {np.min(np.abs(1 - np.exp(eta))):.6f}   (1 - 1/e = {1 - np.exp(-1):.6f})")
    g = np.linspace(1, 6, 400)[:, None]*np.exp(1j*np.linspace(-np.pi, np.pi, 1441))[None, :]
    g = g[np.abs(g.imag) <= np.pi]
    print(f"    min over 1 <= |eta| <= 6, |Im eta| <= pi of |1 - e^eta| = {np.min(np.abs(1 - np.exp(g))):.6f}")
    c = 0.3; m = 1 + np.exp(-c)*np.cos(c); lam = np.exp(c)/m
    print(f"    c = {c}: m_c = 1 + e^(-c) cos c = {m:.4f},  lambda = e^c/m_c = {lam:.4f},  B = 1/(c (1 - lambda)) = {1/(c*(1 - lam)):.2f}")
    c3 = 2*np.exp(-0.5)*np.sin(0.5); L = 2 + 1/c3
    print(f"    proved lower bound c3 = 2 e^(-1/2) sin(1/2) = {c3:.4f};  L = 2 + 1/c3 = {L:.4f};  with c1 = 0.04:  L c1 = {L*0.04:.4f} < kappa - 0.01 = 0.175,  ln 2 - L c1 = {np.log(2) - L*0.04:.4f}")
    print(f"    pi cos(70 deg) - c1 = {np.pi*np.cos(np.radians(70)) - 0.04:.4f}  (distance from the strip around the contour to the points i pi (2m+1); must be >= 1)")
